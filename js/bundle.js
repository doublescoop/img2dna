// img2dna bundle — single non-module script so it works from file:// too.
// All logic in one IIFE. Mirrors img2dna_soul_v11_GR.ipynb exactly.

(() => {
  "use strict";

  // ══════════════════════════════ CODEC ══════════════════════════════════
  const GR_K = 5;
  const MARKER = "TCCG";
  const RESTRICTION_SITES = [
    "GAATTC", "GGTACC", "CCATGG", "GAGCTC", "CGTCTC", "GAGACG",
    "GGTCTC", "GAGACC", "GCGATG", "CATCGC", "GAAGAC", "GTCTTC",
  ];

  function grBitLength(rl, k) { const n = rl - 1; return (n >> k) + 1 + k; }

  function grEncodeRuns(runs, k) {
    const parts = [];
    const mask = (1 << k) - 1;
    for (const rl of runs) {
      const n = rl - 1;
      const q = n >> k;
      const r = n & mask;
      parts.push("1".repeat(q));
      parts.push("0");
      if (k > 0) parts.push(r.toString(2).padStart(k, "0"));
    }
    return parts.join("");
  }

  function grDecodeOne(bits, pos, k) {
    let q = 0; const L = bits.length;
    while (pos < L && bits[pos] === "1") { q++; pos++; }
    if (pos >= L) throw new Error("EOF in unary");
    pos++;
    let r = 0;
    if (k > 0) {
      if (pos + k > L) throw new Error("EOF in remainder");
      r = parseInt(bits.substring(pos, pos + k), 2);
      pos += k;
    }
    return [((q << k) | r) + 1, pos];
  }

  function intToTrits(n) {
    if (n === 0n) return [0];
    const out = [];
    while (n > 0n) { const r = Number(n % 3n); n = n / 3n; out.push(r); }
    return out.reverse();
  }
  function tritsToInt(trits) {
    let n = 0n;
    for (const t of trits) n = n * 3n + BigInt(t);
    return n;
  }
  function bitsToBigInt(bits) { return bits.length === 0 ? 0n : BigInt("0b" + bits); }
  function bigIntToBits(n) { return n === 0n ? "0" : n.toString(2); }

  const FSM = {
    A: { 0: "C", 1: "G", 2: "T" }, C: { 0: "G", 1: "T", 2: "A" },
    G: { 0: "T", 1: "A", 2: "C" }, T: { 0: "A", 1: "C", 2: "G" },
  };
  const FSM_REV = {};
  for (const [p, m] of Object.entries(FSM)) {
    FSM_REV[p] = {};
    for (const [t, c] of Object.entries(m)) FSM_REV[p][c] = Number(t);
  }
  const FIRST_BASE = { 0: "A", 1: "C", 2: "G" };
  const FIRST_BASE_REV = { A: 0, C: 1, G: 2 };

  function tritsToDna(trits) {
    const out = []; let prev = null;
    for (const t of trits) {
      const n = prev === null ? FIRST_BASE[t] : FSM[prev][t];
      out.push(n); prev = n;
    }
    return out.join("");
  }
  function dnaToTrits(dna) {
    const out = []; let prev = null;
    for (const b of dna) {
      let t;
      if (prev === null) { t = FIRST_BASE_REV[b]; if (t === undefined) throw new Error(`bad first base ${b}`); }
      else { t = FSM_REV[prev][b]; if (t === undefined) throw new Error(`bad transition ${prev}->${b}`); }
      out.push(t); prev = b;
    }
    return out;
  }

  function xorWithSalt(bits, salt) {
    if (salt === 0) return bits;
    const key = salt.toString(2).padStart(16, "0");
    const out = new Array(bits.length);
    for (let i = 0; i < bits.length; i++) {
      out[i] = ((bits.charCodeAt(i) - 48) ^ (key.charCodeAt(i % 16) - 48)) ? "1" : "0";
    }
    return out.join("");
  }

  function isDnaSafe(dna) {
    let gc = 0;
    for (const c of dna) if (c === "C" || c === "G") gc++;
    const pct = (gc / dna.length) * 100;
    if (pct < 40 || pct > 60) return false;
    for (const s of RESTRICTION_SITES) if (dna.includes(s)) return false;
    return true;
  }

  function buildChunkBody(chunkRows, k) {
    const parts = [];
    for (const { startBit, runs } of chunkRows) {
      parts.push(startBit.toString());
      parts.push(grEncodeRuns(runs, k));
    }
    return parts.join("");
  }
  function chunkBitsToDna(body, salt) {
    const scrambled = xorWithSalt(body, salt);
    const full = "1" + salt.toString(2).padStart(16, "0") + scrambled;
    return tritsToDna(intToTrits(bitsToBigInt(full)));
  }
  function encodeChunkSearchSalt(chunkRows, k, maxSalts = 65536) {
    const body = buildChunkBody(chunkRows, k);
    for (let salt = 0; salt < maxSalts; salt++) {
      const dna = chunkBitsToDna(body, salt);
      if (isDnaSafe(dna)) return { dna, salt };
    }
    throw new Error("all salts failed");
  }
  function decodeChunkDna(chunkDna, W, k) {
    if (!chunkDna) return [];
    const trits = dnaToTrits(chunkDna);
    const bits = bigIntToBits(tritsToInt(trits));
    const g = bits.indexOf("1");
    if (g < 0 || bits.length < g + 1 + 16) return [];
    const afterGuard = bits.substring(g + 1);
    const salt = parseInt(afterGuard.substring(0, 16), 2);
    const scrambled = afterGuard.substring(16);
    const body = xorWithSalt(scrambled, salt);

    const rows = [];
    let pos = 0;
    while (pos < body.length) {
      const startBit = Number(body[pos]); pos++;
      const runs = []; let total = 0;
      while (total < W) {
        if (pos >= body.length) break;
        try {
          const [v, newPos] = grDecodeOne(body, pos, k);
          if (v > W * 2) break;
          runs.push(v); total += v; pos = newPos;
        } catch (e) { break; }
      }
      rows.push({ startBit, runs });
    }
    return rows;
  }

  function imageToRows(bw, W, H) {
    const rows = [];
    for (let r = 0; r < H; r++) {
      const off = r * W;
      const startBit = bw[off] === 0 ? 0 : 1;
      const runs = [];
      let cur = bw[off], length = 0;
      for (let c = 0; c < W; c++) {
        const p = bw[off + c];
        if (p === cur) length++;
        else { runs.push(length); cur = p; length = 1; }
      }
      runs.push(length);
      rows.push({ startBit, runs });
    }
    return rows;
  }
  function packRowsIntoChunks(rowsData, k, targetBodyBits) {
    const chunks = [];
    let cur = [], bits = 0;
    for (const row of rowsData) {
      const rb = 1 + row.runs.reduce((s, r) => s + grBitLength(r, k), 0);
      if (cur.length && bits + rb > targetBodyBits) { chunks.push(cur); cur = []; bits = 0; }
      cur.push(row); bits += rb;
    }
    if (cur.length) chunks.push(cur);
    return chunks;
  }
  function paintRows(rows, W, H) {
    const canvas = new Uint8Array(H * W).fill(255);
    const n = Math.min(rows.length, H);
    for (let r = 0; r < n; r++) {
      const { startBit, runs } = rows[r];
      let val = startBit === 0 ? 0 : 255, col = 0;
      for (const rl of runs) {
        const eff = Math.min(rl, W - col);
        if (eff <= 0) break;
        canvas.fill(val, r * W + col, r * W + col + eff);
        col += eff; val = 255 - val;
        if (col >= W) break;
      }
    }
    return canvas;
  }
  function decodeDnaSoup(text, W, k, marker = MARKER) {
    const cleaned = text.toUpperCase().replace(/[^ACGT]/g, "");
    const pieces = cleaned.split(marker).filter(Boolean);
    const allRows = [], perChunk = [];
    for (const p of pieces) {
      let rows = [];
      try { rows = decodeChunkDna(p, W, k); } catch { rows = []; }
      perChunk.push(rows.length);
      allRows.push(...rows);
    }
    return { allRows, perChunk };
  }

  // ══════════════════════════════ TABS ═══════════════════════════════════
  const tabs = document.querySelectorAll(".tab");
  const panels = document.querySelectorAll(".panel");

  function activate(name) {
    tabs.forEach((t) => t.classList.toggle("active", t.dataset.tab === name));
    panels.forEach((p) => p.classList.toggle("active", p.dataset.panel === name));
    history.replaceState(null, "", `#${name}`);
  }
  tabs.forEach((t) => t.addEventListener("click", () => activate(t.dataset.tab)));
  const initial = location.hash.replace("#", "") || "decode";
  activate(initial === "encode" ? "encode" : "decode");

  // ══════════════════════════════ DECODER UI ═════════════════════════════
  const decIn = document.getElementById("dec-input");
  const decWidth = document.getElementById("dec-width");
  const decCanvas = document.getElementById("dec-canvas");
  const decStats = document.getElementById("dec-stats");
  const decFile = document.getElementById("dec-file");
  const decClear = document.getElementById("dec-clear");
  let lastDecKey = "";

  function renderDecode() {
    const text = decIn.value;
    const W = parseInt(decWidth.value, 10) || 200;
    const key = `${W}::${text}`;
    if (key === lastDecKey) return;
    lastDecKey = key;

    if (!text.trim()) {
      decCanvas.width = W; decCanvas.height = 1;
      const ctx = decCanvas.getContext("2d");
      ctx.fillStyle = "white"; ctx.fillRect(0, 0, W, 1);
      decStats.textContent = `Paste DNA to render. Expected marker: ${MARKER}. Expected width: ${W}.`;
      return;
    }
    const t0 = performance.now();
    const { allRows, perChunk } = decodeDnaSoup(text, W, GR_K);
    const elapsed = (performance.now() - t0).toFixed(0);

    if (allRows.length === 0) {
      decCanvas.width = W; decCanvas.height = 1;
      const ctx = decCanvas.getContext("2d");
      ctx.fillStyle = "white"; ctx.fillRect(0, 0, W, 1);
      decStats.textContent =
        `No rows decoded from ${perChunk.length} candidate chunk(s). ` +
        `Check W=${W}, GR_K=${GR_K}, MARKER=${MARKER}.`;
      return;
    }
    const H = allRows.length;
    const bw = paintRows(allRows, W, H);
    decCanvas.width = W; decCanvas.height = H;
    const ctx = decCanvas.getContext("2d");
    const imgData = ctx.createImageData(W, H);
    const d = imgData.data;
    for (let i = 0, p = 0; p < bw.length; p++, i += 4) {
      const v = bw[p]; d[i] = v; d[i + 1] = v; d[i + 2] = v; d[i + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
    decStats.textContent =
      `${perChunk.length} chunk(s) found  ·  rows per chunk: [${perChunk.join(", ")}]  ·  ` +
      `${H} rows painted × ${W} cols  ·  ${elapsed} ms`;
  }
  decIn.addEventListener("input", renderDecode);
  decWidth.addEventListener("input", renderDecode);
  decFile.addEventListener("change", async (e) => {
    const file = e.target.files[0]; if (!file) return;
    decIn.value = await file.text();
    renderDecode();
  });
  decClear.addEventListener("click", () => { decIn.value = ""; renderDecode(); });
  renderDecode();

  // ══════════════════════════════ ENCODER UI ═════════════════════════════
  const encFile = document.getElementById("enc-file");
  const encContrast = document.getElementById("enc-contrast");
  const encThreshold = document.getElementById("enc-threshold");
  const encMinWidth = document.getElementById("enc-min-width");
  const encTargetBits = document.getElementById("enc-target-bits");
  const encContrastVal = document.getElementById("enc-contrast-val");
  const encThresholdVal = document.getElementById("enc-threshold-val");
  const encMinWidthVal = document.getElementById("enc-min-width-val");
  const encTargetBitsVal = document.getElementById("enc-target-bits-val");
  const encPreview = document.getElementById("enc-preview");
  const encBtn = document.getElementById("enc-encode-btn");
  const encStatus = document.getElementById("enc-status");
  const encResults = document.getElementById("enc-results");
  const encSoup = document.getElementById("enc-soup");
  const encFasta = document.getElementById("enc-fasta");
  const encDlSoup = document.getElementById("enc-download-soup");
  const encDlFasta = document.getElementById("enc-download-fasta");
  const encCopySoup = document.getElementById("enc-copy-soup");

  let loadedImage = null, bwProcessed = null, encH = 0, encW = 0;

  function syncLabels() {
    encContrastVal.textContent = encContrast.value;
    encThresholdVal.textContent = encThreshold.value;
    encMinWidthVal.textContent = encMinWidth.value;
    encTargetBitsVal.textContent = encTargetBits.value;
  }
  [encContrast, encThreshold, encMinWidth, encTargetBits].forEach((el) =>
    el.addEventListener("input", () => { syncLabels(); if (loadedImage) renderEncPreview(); }));
  syncLabels();

  encFile.addEventListener("change", async (e) => {
    const file = e.target.files[0]; if (!file) return;
    const url = URL.createObjectURL(file);
    const img = new Image();
    img.onload = () => { loadedImage = img; renderEncPreview(); URL.revokeObjectURL(url); };
    img.src = url;
  });

  function renderEncPreview() {
    if (!loadedImage) return;
    const minW = parseInt(encMinWidth.value, 10);
    const contrast = parseFloat(encContrast.value);
    const threshold = parseInt(encThreshold.value, 10);
    const scale = minW / loadedImage.naturalWidth;
    const tw = minW;
    const th = Math.max(1, Math.round(loadedImage.naturalHeight * scale));
    encPreview.width = tw; encPreview.height = th;
    const ctx = encPreview.getContext("2d", { willReadFrequently: true });
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(loadedImage, 0, 0, tw, th);
    const imgData = ctx.getImageData(0, 0, tw, th);
    const d = imgData.data;
    const midgray = 128;
    const bw = new Uint8Array(tw * th);
    for (let i = 0, p = 0; i < d.length; i += 4, p++) {
      let g = 0.2989 * d[i] + 0.587 * d[i + 1] + 0.114 * d[i + 2];
      g = (g - midgray) * contrast + midgray;
      g = Math.max(0, Math.min(255, g));
      const v = g < threshold ? 0 : 255;
      bw[p] = v;
      d[i] = d[i + 1] = d[i + 2] = v; d[i + 3] = 255;
    }
    ctx.putImageData(imgData, 0, 0);
    bwProcessed = bw; encW = tw; encH = th;
    encBtn.disabled = false;
    encStatus.textContent = `Preview: ${encH} rows × ${encW} cols (${(encH * encW).toLocaleString()} pixels).`;
  }

  encBtn.addEventListener("click", async () => {
    if (!bwProcessed) return;
    encStatus.textContent = "Encoding…";
    encBtn.disabled = true; encResults.hidden = true;
    await new Promise((r) => setTimeout(r, 10));
    const targetBits = parseInt(encTargetBits.value, 10);
    const t0 = performance.now();
    const rows = imageToRows(bwProcessed, encW, encH);
    const chunks = packRowsIntoChunks(rows, GR_K, targetBits);
    const dnaChunks = [], salts = [], rowsPerChunk = [];
    try {
      for (let i = 0; i < chunks.length; i++) {
        const { dna, salt } = encodeChunkSearchSalt(chunks[i], GR_K);
        dnaChunks.push(dna); salts.push(salt); rowsPerChunk.push(chunks[i].length);
        encStatus.textContent = `Chunk ${i + 1}/${chunks.length} (salt=${salt})…`;
        await new Promise((r) => setTimeout(r, 0));
      }
    } catch (e) {
      encStatus.textContent = `Error: ${e.message}. Try lowering target body bits.`;
      encBtn.disabled = false; return;
    }
    const elapsed = ((performance.now() - t0) / 1000).toFixed(2);
    const soup = dnaChunks.map((d) => MARKER + d).join("");
    const header = `# W=${encW} GR_K=${GR_K} MARKER=${MARKER} H=${encH} chunks=${dnaChunks.length}`;
    const fasta = header + "\n" +
      dnaChunks.map((d, i) => `>chunk_${i} salt=${salts[i]} rows=${rowsPerChunk[i]}\n${MARKER}${d}`).join("\n") + "\n";
    encSoup.value = header + "\n" + soup;
    encFasta.value = fasta;
    encResults.hidden = false; encBtn.disabled = false;
    encStatus.textContent = `Done. ${soup.length.toLocaleString()} DNA bases in ${dnaChunks.length} chunks (${elapsed}s). ` +
      `${(encH * encW / soup.length).toFixed(2)} pixels/base.`;
  });

  function downloadText(filename, text) {
    const blob = new Blob([text], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename;
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  }
  encDlSoup.addEventListener("click", () => downloadText("dna_v11_soup.txt", encSoup.value));
  encDlFasta.addEventListener("click", () => downloadText("dna_v11_chunks.fasta", encFasta.value));
  encCopySoup.addEventListener("click", async () => {
    await navigator.clipboard.writeText(encSoup.value);
    encCopySoup.textContent = "Copied!";
    setTimeout(() => { encCopySoup.textContent = "Copy soup"; }, 1200);
  });
})();

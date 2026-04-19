# img2dna

Static web app for encoding a black-and-white image as DNA (and decoding DNA back to an image). Mirrors the v11 notebook exactly — same constants, same format, roundtrips.

## Host on GitHub Pages

1. Create a new GitHub repo (e.g. `img2dna`) and push the contents of this folder:
   ```
   git init
   git add .
   git commit -m "initial"
   git branch -M main
   git remote add origin git@github.com:<you>/img2dna.git
   git push -u origin main
   ```
2. In the repo on github.com → **Settings → Pages**:
   - **Source**: Deploy from a branch
   - **Branch**: `main` / `/ (root)`
3. Wait ~30s. Site lives at `https://<you>.github.io/img2dna/`.

No build step. No dependencies. The site is pure HTML + ES modules.

## Run locally

Any static server. Examples:

```bash
python3 -m http.server 8000   # then open http://localhost:8000
```
or
```bash
npx serve .
```

## Format

| Constant | Value |
|---|---|
| `W` | image width in pixels (set at encode time) |
| `GR_K` | `5` |
| `MARKER` | `TCCG` |

Per strand: `TCCG` + trit-encoded (guard=1b, salt=16b, XOR'd GR body).
Body = interleaved `(start_bit, GR runs summing to W)` per row.

## Decoder spec (for biologists)

1. Keep only `A/C/G/T` (uppercase).
2. Split on `TCCG`.
3. For each piece:
   - Convert DNA → trits using the FSM.
   - Trits → bits.
   - Skip leading zeros; first `1` is the guard.
   - Next 16 bits = salt.
   - XOR remaining bits with salt (repeating every 16 bits).
   - Loop: read 1 bit (row start color), GR-decode runs until they sum to `W` → paint row.
4. Paint rows top-to-bottom on a white canvas.

Partial DNA gives partial picture — missing chunks leave white rows; unreadable chunks are skipped.

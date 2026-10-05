# Vendored browser libraries

Every file in this directory was fetched over HTTPS from `https://unpkg.com/`
with `curl`, on 2026-10-05, at the exact version recorded below. No file here
was hand-written, minified further, or re-fetched at a floating version.

**Total: 109,383 bytes — 103,535 of JS plus 5,848 of licence text — which is
106.82 KiB, 0.1043 MiB, or 1.74 % of the 6 MB budget.**

`du -sk` on this directory reports 148 KiB, because every 731-byte licence file
occupies a 4 KiB block. The sum of the file sizes is the number to quote, and it
is the sum of the sizes.

| file | version | bytes | license | SHA-256 |
|---|---|---|---|---|
| `d3-array.min.js` | 3.2.4 | 17,204 | ISC | `80aa70d0cd17dabddf6d056494ea17926a45a69da8b7850220aace331bad671d` |
| `d3-axis.min.js` | 3.0.0 | 3,118 | ISC | `e3748d1f430223fbdc207d67e0d8962031e8920f4a7c6cb572241a30b1e08dc7` |
| `d3-format.min.js` | 3.1.0 | 5,218 | ISC | `7d053f71a135100128802b6010ab58d1351ca412e2d7846c2f8c6fe155a66370` |
| `d3-interpolate.min.js` | 3.0.1 | 7,863 | ISC | `bfc321e4c3f3b3aadc88cfe15ccb5e443abfeadef8b75c65b41c33a4d78a98ae` |
| `d3-scale.min.js` | 4.0.2 | 15,728 | ISC | `e76a84839ffba3b94fef22ea1e39da8398fa0e039c7e6a0a93b7d938dbd50632` |
| `d3-selection.min.js` | 3.0.0 | 13,522 | ISC | `45daab9cf677901bcae102f3f23ca2930db3c0fb8ff9e3dbed087d9c4de921ca` |
| `d3-shape.min.js` | 3.2.0 | 30,898 | ISC | `49236896e989f251dd36679ea4a879538e9ca2ece3a90a82740f953cac4a3fa5` |
| `d3-zoom.min.js` | 3.0.0 | 9,984 | ISC | `4fdce9b830b78225c58e27258296084253c41a5eaf8cad9bd9ccfe3025daf1cf` |

The eight `*.LICENSE` files beside them are the upstream `LICENSE` files for the
same eight versions, 731 bytes each (5,848 bytes total), also fetched from
unpkg and byte-identical to upstream.

## Why this set and not a bundle

Each file is a UMD build, so the browser loads them as plain `<script>` tags with
no bundler, no `node_modules`, and no build step. They were chosen one at a time
against the four visualisations V1–V4 actually need:

| needed for | module |
|---|---|
| every SVG/DOM operation in all four | `d3-selection` |
| pan + zoom on the tree and the matrix | `d3-zoom` |
| linear/log/band scales: Manhattan x-axis, QQ y-axis, matrix colour, tree x | `d3-scale` |
| tick rendering for the Manhattan and QQ axes | `d3-axis` |
| the radial/dendrogram link generator for the tree | `d3-shape` |
| number formatting for p-values and identity percentages | `d3-format` |
| transitive: `scale` needs `array` and `interpolate` | `d3-array`, `d3-interpolate` |

Deliberately **not** vendored:

- **A table library.** The 900-isolate grid is a custom virtualised canvas
  component. A DOM table library cannot hold 14,400 cells at 60fps, and no
  general-purpose grid handles the "not assessed" badge vocabulary that V2 needs.
  See `DESIGN.md` § D7.
- **`d3-force` / `d3-delaunay`.** The co-occurrence view (V4) draws a fixed
  circular layout. A force simulation on a 900-isolate cohort would settle
  differently on every render, which makes the picture non-reproducible — a
  diagram that moves between two views of the same run is a diagram that cannot
  be cited.
- **`d3-hierarchy`.** The tree endpoint (V3) already returns structural JSON with
  server-assigned node ids; a client-side layout pass would have to re-derive the
  same ids the server assigned, which is precisely the D6 failure the design
  forbids.

## Verification performed

```
curl -sS -m 30 -w "%{http_code}" -o <file> https://unpkg.com/<module>@<version>/dist/<file>.min.js
# -> HTTP 200 for all eight
node --check <file>     # -> exit 0 for all eight
shasum -a 256 <file>    # -> the digests in the table above
```

Each downloaded file begins with a banner naming its own version
(`// https://d3js.org/d3-selection/ v3.0.0 Copyright 2010-2021 Mike Bostock`),
which is how the bytes were checked to be the release claimed and not an error
page.

## Licence

All eight modules are © Mike Bostock and released under the ISC licence, which
permits redistribution. The licence text is reproduced in full in each
`*.LICENSE` file. Nothing here is redistributed under any other term.
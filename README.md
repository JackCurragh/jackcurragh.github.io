# CitationLinker

CitationLinker is a fully client-side web app for converting Paperpile-generated external DOCX citation links into internal DOCX navigation links.

## What it does

- Accepts a Paperpile DOCX by drag-and-drop or file picker.
- Detects Paperpile citation and bibliography hyperlinks in the DOCX XML locally in the browser.
- Maps Paperpile `/c/<doc>/<ids>` citations to matching `/b/<doc>/<id>` bibliography entries.
- Rewrites bibliography hyperlinks as DOCX bookmarks and citation hyperlinks as internal DOCX anchors.
- Preserves unrelated document content and hyperlinks.
- Downloads a modified DOCX without uploading the manuscript anywhere.

Semicolon-separated multi-reference citation hyperlinks are split into separate DOCX hyperlink runs when the text and Paperpile ID counts agree. Other multi-reference groups are reported for review rather than guessed.

## Privacy

All processing happens in the browser. No manuscript content is uploaded or transmitted.

## Development

```bash
npm install
npm run dev
```

## GitHub Pages

This repo is deployed as a static site. GitHub Pages should publish the `dist/` output produced by `npm run build`.

The source `index.html` is a Vite entry point, so opening it directly via `file://` will fail. Use `npm run dev` locally, or view the deployed GitHub Pages site.

## Build

```bash
npm run build
```

## Tests

```bash
npm test
```

## Notes

- Multi-reference citation clusters are handled conservatively.
- Ambiguous clusters are reported rather than guessed.
- The app is designed for GitHub Pages deployment as a static site.

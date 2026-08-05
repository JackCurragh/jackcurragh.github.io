# Jack Tierney - Personal Website

This repository contains Jack Tierney's personal website and a collection of standalone browser apps for research and translation analysis.

## Site Structure

```text
├── index.html              # Home page
├── research.html           # Research and projects
├── publications.html       # Publications
├── tools.html              # Tools
├── apps/                   # Interactive apps and visualizations
│   ├── index.html          # App gallery
│   └── citation-linker/    # Vite-built Paperpile DOCX linker
└── assets/                 # Shared site assets
```

## CitationLinker

CitationLinker is a fully client-side app at `apps/citation-linker/`. It converts Paperpile-generated external DOCX citation links into internal DOCX navigation links without uploading the manuscript.

```bash
npm install
npm test
npm run build
```

The Vite build writes the app to `dist/apps/citation-linker/`. The GitHub Pages workflow assembles the existing static website alongside it.

## GitHub Pages

Set **Pages -> Build and deployment -> Source** to **GitHub Actions**. The deployment workflow publishes the existing website at `/` and the apps gallery at `/apps/`.

## License

Content: Copyright 2024 Jack Tierney
Code: MIT License (see LICENSE file)

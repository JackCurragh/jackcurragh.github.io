# Jack Tierney - Personal Website

A minimal personal website for showcasing research, publications, and interactive bioinformatics visualizations.

## Site Structure

```
├── index.html          # Home page with bio and overview
├── research.html       # Publications and ongoing projects
├── demos/              # Interactive visualizations
│   ├── index.html      # Gallery of demos
│   └── example-demo/   # Template for new demos
├── assets/
│   ├── css/
│   │   └── style.css   # Main stylesheet
│   └── img/
│       └── jack.jpeg   # Profile photo
```

## Adding New Demos

Each demo should be a standalone HTML page in its own directory under `demos/`:

1. Create a new directory: `demos/your-demo-name/`
2. Add an `index.html` file with your visualization
3. Include any necessary JavaScript libraries (D3.js, Three.js, etc.)
4. Link to the demo from `demos/index.html`

### Example Demo Structure

```
demos/
└── ribosome-visualization/
    ├── index.html
    ├── script.js
    ├── style.css
    └── data/
        └── sample-data.json
```

## Technologies

- Pure HTML/CSS/JavaScript (no build system required)
- Hosted on GitHub Pages
- Designed for adding interactive visualizations using:
  - D3.js for data visualization
  - Three.js for 3D graphics
  - Canvas API for custom animations

## Development

Simply edit the HTML/CSS files and push to GitHub. The site will automatically update via GitHub Pages.

## License

Content: © 2024 Jack Tierney
Code: MIT License (see LICENSE file)

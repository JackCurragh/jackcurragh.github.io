// Theme Switcher
// Manages theme selection and persistence across pages

(function() {
    'use strict';

    const THEMES = ['academic', 'terminal', 'dataviz', 'brutalist', 'retro'];
    const DEFAULT_THEME = 'academic';
    const STORAGE_KEY = 'preferred-theme';

    // Get current theme from localStorage or default
    function getCurrentTheme() {
        const stored = localStorage.getItem(STORAGE_KEY);
        return THEMES.includes(stored) ? stored : DEFAULT_THEME;
    }

    // Apply theme to body element
    function applyTheme(theme) {
        // Remove all theme classes
        THEMES.forEach(t => {
            document.body.classList.remove(`theme-${t}`);
        });

        // Add selected theme class
        document.body.classList.add(`theme-${theme}`);

        // Store preference
        localStorage.setItem(STORAGE_KEY, theme);

        // Update selector if it exists
        const selector = document.getElementById('theme-selector');
        if (selector) {
            selector.value = theme;
        }
    }

    // Create theme selector dropdown
    function createThemeSelector() {
        const container = document.createElement('div');
        container.className = 'theme-switcher';
        container.style.cssText = `
            position: fixed;
            top: 1rem;
            right: 1rem;
            z-index: 1000;
            padding: 0.75rem 1rem;
            background: rgba(255, 255, 255, 0.95);
            border: 2px solid #000;
            border-radius: 4px;
            box-shadow: 3px 3px 0 rgba(0, 0, 0, 0.2);
            font-family: inherit;
        `;

        const label = document.createElement('label');
        label.textContent = 'STYLE: ';
        label.style.cssText = `
            font-size: 0.85rem;
            font-weight: bold;
            margin-right: 0.5rem;
            letter-spacing: 0.05em;
        `;

        const select = document.createElement('select');
        select.id = 'theme-selector';
        select.style.cssText = `
            font-size: 1rem;
            padding: 0.4rem 0.75rem;
            border: 2px solid #000;
            background: #fff;
            font-weight: 600;
            cursor: pointer;
            font-family: inherit;
        `;

        THEMES.forEach(theme => {
            const option = document.createElement('option');
            option.value = theme;
            option.textContent = theme.charAt(0).toUpperCase() + theme.slice(1);
            select.appendChild(option);
        });

        select.value = getCurrentTheme();

        select.addEventListener('change', (e) => {
            applyTheme(e.target.value);
        });

        container.appendChild(label);
        container.appendChild(select);

        return container;
    }

    // Initialize on page load
    function init() {
        // Apply saved theme immediately
        const currentTheme = getCurrentTheme();
        applyTheme(currentTheme);

        // Add theme selector to page
        document.addEventListener('DOMContentLoaded', () => {
            const selector = createThemeSelector();
            document.body.appendChild(selector);
        });
    }

    init();
})();

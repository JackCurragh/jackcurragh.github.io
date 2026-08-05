import { defineConfig } from "vite";

export default defineConfig({
  base: "./",
  build: {
    outDir: "dist",
    emptyOutDir: false,
    assetsDir: "apps/citation-linker/assets",
    rollupOptions: {
      input: "apps/citation-linker/index.html",
    },
  },
  test: {
    environment: "jsdom",
  },
});

import { cpSync, mkdirSync, rmSync } from "node:fs";

rmSync("dist", { force: true, recursive: true });
mkdirSync("dist", { recursive: true });

for (const file of ["index.html", "research.html", "tools.html", "publications.html"]) {
  cpSync(file, `dist/${file}`);
}

cpSync("assets", "dist/assets", { recursive: true });
cpSync("apps", "dist/apps", { recursive: true });

import type { AppState } from "./state";
import { summarizeProcessing } from "./validation";
import type { CitationProblem } from "../model/types";

export interface AppActions {
  focusResolver(problemKey: string): void;
}

export interface AppUi {
  root: HTMLElement;
  statusEl: HTMLElement;
  reportEl: HTMLElement;
  downloadEl: HTMLAnchorElement;
  problemsDownloadEl: HTMLAnchorElement;
  problemsSummaryEl: HTMLElement;
  problemsListEl: HTMLOListElement;
  resolverPanelEl: HTMLElement;
  resolverHeaderEl: HTMLElement;
  resolverMetaEl: HTMLElement;
  resolverFragmentsEl: HTMLElement;
  resolverEmptyEl: HTMLElement;
  errorEl: HTMLElement;
  fileLabelEl: HTMLElement;
}

export function createAppShell(): AppUi {
  const root = document.createElement("div");
  root.className = "app-shell";
  root.innerHTML = `
    <section class="hero">
      <div class="eyebrow">Offline Paperpile tooling</div>
      <h1>CitationLinker</h1>
      <p class="lede">Convert Paperpile-generated external citation links into internal DOCX navigation entirely in your browser.</p>
      <p class="privacy">Your manuscript is processed entirely in your browser. It is not uploaded or transmitted anywhere.</p>
    </section>
    <section class="panel dropzone" id="dropzone">
      <input id="file-input" type="file" accept="application/vnd.openxmlformats-officedocument.wordprocessingml.document,.docx" hidden />
      <div class="dropzone-copy">
        <strong>Drop your Paperpile DOCX here</strong>
        <span>DOCX keeps the citation and bibliography link boundaries intact.</span>
      </div>
      <button class="primary" id="choose-button" type="button">Choose DOCX</button>
    </section>
    <section class="panel controls">
      <div class="row">
        <div>
          <div class="label">Selected file</div>
          <div id="file-label" class="value">None</div>
        </div>
        <div class="actions">
          <button class="secondary" id="process-button" type="button" disabled>Analyze and rewrite</button>
          <a id="download-link" class="download" download="citation-linker-output.docx" hidden>Download DOCX</a>
          <a id="problems-download-link" class="download" download="citation-linker-problems.json" hidden>Download problems</a>
        </div>
      </div>
      <div id="status" class="status">Drop a Paperpile DOCX to begin.</div>
      <div id="error" class="error" hidden></div>
    </section>
    <section class="panel report">
      <div class="label">Validation report</div>
      <pre id="report">No Paperpile DOCX processed yet.</pre>
    </section>
    <section class="panel report">
      <div class="label">Resolver</div>
      <div id="resolver-empty" class="status">No unresolved citation clusters.</div>
      <div id="resolver-panel" class="resolver-panel" hidden>
        <div class="resolver-help">Paperpile IDs are shown beside the citation text extracted from the DOCX. Unresolved groups were not rewritten.</div>
        <div id="resolver-header" class="resolver-header"></div>
        <div id="resolver-meta" class="resolver-meta"></div>
        <div id="resolver-fragments" class="resolver-fragments"></div>
      </div>
      <ol id="problems-list" class="problem-list"></ol>
    </section>
  `;

  document.body.appendChild(root);

  return {
    root,
    statusEl: root.querySelector("#status") as HTMLElement,
    reportEl: root.querySelector("#report") as HTMLElement,
    downloadEl: root.querySelector("#download-link") as HTMLAnchorElement,
    problemsDownloadEl: root.querySelector("#problems-download-link") as HTMLAnchorElement,
    problemsSummaryEl: root.querySelector("#resolver-empty") as HTMLElement,
    problemsListEl: root.querySelector("#problems-list") as HTMLOListElement,
    resolverPanelEl: root.querySelector("#resolver-panel") as HTMLElement,
    resolverHeaderEl: root.querySelector("#resolver-header") as HTMLElement,
    resolverMetaEl: root.querySelector("#resolver-meta") as HTMLElement,
    resolverFragmentsEl: root.querySelector("#resolver-fragments") as HTMLElement,
    resolverEmptyEl: root.querySelector("#resolver-empty") as HTMLElement,
    errorEl: root.querySelector("#error") as HTMLElement,
    fileLabelEl: root.querySelector("#file-label") as HTMLElement,
  };
}

export function applyState(ui: AppUi, state: AppState, actions: AppActions): void {
  ui.statusEl.textContent = state.status;
  ui.fileLabelEl.textContent = state.fileName ?? "None";
  ui.reportEl.textContent = state.report ? summarizeProcessing(state.report) : "No Paperpile DOCX processed yet.";
  renderProblems(ui, state, actions);

  if (state.error) {
    ui.errorEl.hidden = false;
    ui.errorEl.textContent = state.error;
  } else {
    ui.errorEl.hidden = true;
    ui.errorEl.textContent = "";
  }

  if (state.outputUrl) {
    ui.downloadEl.hidden = false;
    ui.downloadEl.href = state.outputUrl;
  } else {
    ui.downloadEl.hidden = true;
    ui.downloadEl.removeAttribute("href");
  }

  if (state.problemsUrl) {
    ui.problemsDownloadEl.hidden = false;
    ui.problemsDownloadEl.href = state.problemsUrl;
  } else {
    ui.problemsDownloadEl.hidden = true;
    ui.problemsDownloadEl.removeAttribute("href");
  }

  const hasProblem = Boolean(state.report?.problems?.length);
  ui.resolverEmptyEl.hidden = hasProblem;
  ui.resolverPanelEl.hidden = !state.resolverFocusKey;

  const focusedProblem = state.report?.problems.find((problem) => problem.problemKey === state.resolverFocusKey);
  renderResolver(ui, state, focusedProblem, actions);
}

function renderProblems(ui: AppUi, state: AppState, actions: AppActions): void {
  const problems = state.report?.problems ?? [];
  if (problems.length === 0) {
    ui.problemsSummaryEl.textContent = "No unresolved citation clusters.";
    ui.problemsListEl.replaceChildren();
    return;
  }

  const reasonCounts = new Map<string, number>();
  for (const problem of problems) {
    reasonCounts.set(problem.code, (reasonCounts.get(problem.code) ?? 0) + 1);
  }

  ui.problemsSummaryEl.textContent = `${problems.length} unresolved citation clusters. ` +
    [...reasonCounts.entries()].map(([reason, count]) => `${describeProblemCode(reason)}: ${count}`).join(", ");

  const items = problems.map((problem) => {
    const item = document.createElement("li");
    item.className = "problem-item";

    const header = document.createElement("div");
    header.className = "problem-header";
    header.textContent = `${problem.pageIndex < 0 ? "DOCX" : `Page ${problem.pageIndex + 1}`} · ${describeProblemCode(problem.code)}`;

    const subheader = document.createElement("div");
    subheader.className = "problem-subheader";
    subheader.textContent = "This Paperpile citation group was not rewritten.";

    if (problem.previewDataUrl) {
      const img = document.createElement("img");
      img.className = "problem-preview";
      img.src = problem.previewDataUrl;
      img.alt = `Page ${problem.pageIndex + 1} citation preview`;
      item.append(img);
    }

    const refs = document.createElement("div");
    refs.className = "problem-refs";
    refs.textContent = `Refs: ${problem.referenceIds.join(", ")}`;

    const fragmentCount = document.createElement("div");
    fragmentCount.className = "problem-refs";
    fragmentCount.textContent = `Fragments: ${problem.fragments.length}`;

    const pairingPreview = document.createElement("div");
    pairingPreview.className = "problem-pairing-preview";
    pairingPreview.textContent = renderPairingPreview(problem);

    const uri = document.createElement("div");
    uri.className = "problem-uri";
    uri.textContent = problem.uri;

    const message = document.createElement("div");
    message.className = "problem-message";
    message.textContent = problem.message;

    const evidence = document.createElement("div");
    evidence.className = "problem-evidence";
    evidence.textContent = problem.evidenceText ? `Extracted text: ${problem.evidenceText}` : "Extracted text unavailable.";

    const focusButton = document.createElement("button");
    focusButton.className = "secondary";
    focusButton.type = "button";
    focusButton.textContent = "Inspect pair";
    focusButton.addEventListener("click", () => actions.focusResolver(problem.problemKey));

    item.append(header, subheader, refs, fragmentCount, pairingPreview, uri, message, evidence, focusButton);
    return item;
  });

  ui.problemsListEl.replaceChildren(...items);
}

function renderResolver(ui: AppUi, state: AppState, problem: CitationProblem | undefined, actions: AppActions): void {
  if (!problem) {
    ui.resolverPanelEl.hidden = true;
    ui.resolverHeaderEl.textContent = "";
    ui.resolverMetaEl.textContent = "";
    ui.resolverFragmentsEl.replaceChildren();
    return;
  }

  ui.resolverPanelEl.hidden = false;
  ui.resolverHeaderEl.textContent = `${problem.pageIndex < 0 ? "DOCX" : `Page ${problem.pageIndex + 1}`} · ${describeProblemCode(problem.code)}`;
  ui.resolverMetaEl.textContent = `${problem.referenceIds.length} Paperpile IDs, ${problem.fragments.length} extracted citation segments`;

  const fragments = problem.fragments.map((fragment) => {
    const row = document.createElement("div");
    row.className = "pairing-row";

    const fragmentLabel = document.createElement("div");
    fragmentLabel.className = "pairing-fragment";
    fragmentLabel.textContent = fragment.text || `Fragment ${fragment.index + 1}`;

    const arrow = document.createElement("div");
    arrow.className = "pairing-arrow";
    arrow.textContent = "↔";

    const reference = document.createElement("div");
    reference.className = "pairing-reference";
    reference.textContent = problem.referenceIds[Math.min(fragment.index, problem.referenceIds.length - 1)] ?? "Unassigned";

    row.append(fragmentLabel, arrow, reference);
    return row;
  });
  ui.resolverFragmentsEl.replaceChildren(...fragments);

}

function renderPairingPreview(problem: CitationProblem): string {
  if (problem.fragments.length === 0 || problem.referenceIds.length === 0) {
    return "No pairing preview available.";
  }

  const pairs = problem.fragments.map((fragment, index) => {
    const refId = problem.referenceIds[Math.min(index, problem.referenceIds.length - 1)];
    return `${fragment.text || `Fragment ${index + 1}`} ↔ ${refId}`;
  });

  return `Suggested pairs: ${pairs.join(" | ")}`;
}

function describeProblemCode(code: string): string {
  switch (code) {
    case "missing-bibliography":
      return "Missing full reference";
    case "no-text-spans":
      return "Could not read the citation text";
    case "docx-multi-reference":
      return "DOCX citation contains multiple references";
    case "segment-count-mismatch":
      return "Citation text does not split cleanly";
    case "text-match-failed":
      return "Could not align citation text";
    case "span-match-failed":
      return "Could not map text spans";
    case "rect-union-failed":
      return "Could not build a clickable area";
    case "no-manual-fragments":
      return "No visible citation fragments";
    case "manual-assignment-incomplete":
      return "Manual assignment not finished";
    case "manual-reference-unassigned":
      return "One fragment is unassigned";
    default:
      return code;
  }
}

import "./style.css";
import { createAppShell, applyState, type AppActions } from "./app/ui";
import { createInitialState, type AppState } from "./app/state";
import { formatError } from "./app/errors";
import { processDocx } from "./docx/process";

const ui = createAppShell();
let selectedFile: File | undefined;
let currentOutputUrl: string | undefined;
let currentProblemsUrl: string | undefined;
let state: AppState = createInitialState();

const dropzone = ui.root.querySelector("#dropzone") as HTMLElement;
const fileInput = ui.root.querySelector("#file-input") as HTMLInputElement;
const chooseButton = ui.root.querySelector("#choose-button") as HTMLButtonElement;
const processButton = ui.root.querySelector("#process-button") as HTMLButtonElement;

function setState(patch: Partial<AppState>): void {
  state = { ...state, ...patch };
  applyState(ui, state, actions);
}

async function runRewrite(): Promise<void> {
  if (!selectedFile) return;

  try {
    if (currentOutputUrl) {
      URL.revokeObjectURL(currentOutputUrl);
      currentOutputUrl = undefined;
    }
    if (currentProblemsUrl) {
      URL.revokeObjectURL(currentProblemsUrl);
      currentProblemsUrl = undefined;
    }

    setState({ status: "Processing locally…", error: undefined, problemsUrl: undefined });
    const { report, outputBytes } = await processDocx(selectedFile, (message) => setState({ status: message }));
    const blob = new Blob([outputBytes], { type: "application/vnd.openxmlformats-officedocument.wordprocessingml.document" });
    currentOutputUrl = URL.createObjectURL(blob);
    const problemsForDownload = report.problems.map(({ previewDataUrl, ...problem }) => problem);
    const problemsBlob = new Blob([JSON.stringify(problemsForDownload, null, 2)], { type: "application/json" });
    currentProblemsUrl = URL.createObjectURL(problemsBlob);
    ui.downloadEl.download = relinkedFileName(selectedFile.name);
    ui.problemsDownloadEl.download = problemsFileName(selectedFile.name);
    setState({
      status: "Finished. Download the rewritten DOCX.",
      report,
      outputUrl: currentOutputUrl,
      problemsUrl: currentProblemsUrl,
      error: undefined,
    });
  } catch (error) {
    setState({
      status: "Processing failed.",
      error: formatError(error),
      outputUrl: undefined,
      problemsUrl: undefined,
    });
  }
}

const actions: AppActions = {
  focusResolver(problemKey) {
    setState({
      resolverFocusKey: problemKey,
    });
  },
};

chooseButton.addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", () => {
  const file = fileInput.files?.[0];
  if (file) {
    selectFile(file);
  }
});

dropzone.addEventListener("dragover", (event) => {
  event.preventDefault();
  dropzone.classList.add("dragging");
});

dropzone.addEventListener("dragleave", () => {
  dropzone.classList.remove("dragging");
});

dropzone.addEventListener("drop", (event) => {
  event.preventDefault();
  dropzone.classList.remove("dragging");
  const file = event.dataTransfer?.files?.[0];
  if (file) {
    selectFile(file);
  }
});

processButton.addEventListener("click", () => {
  void runRewrite();
});

function selectFile(file: File): void {
  selectedFile = file;
  processButton.disabled = false;
  ui.downloadEl.download = relinkedFileName(file.name);
  ui.problemsDownloadEl.download = problemsFileName(file.name);
  setState({
    fileName: file.name,
    status: "Ready to analyze.",
    error: undefined,
    report: undefined,
    outputUrl: undefined,
    problemsUrl: undefined,
    resolverFocusKey: undefined,
  });
}

function relinkedFileName(fileName: string): string {
  return `${fileName.replace(/\.docx$/i, "")}_relinked.docx`;
}

function problemsFileName(fileName: string): string {
  return `${fileName.replace(/\.docx$/i, "")}_problems.json`;
}

applyState(ui, state, actions);

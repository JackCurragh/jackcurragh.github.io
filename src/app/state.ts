import type { ProcessingReport } from "../model/types";

export interface AppState {
  status: string;
  fileName?: string;
  error?: string;
  report?: ProcessingReport;
  outputUrl?: string;
  problemsUrl?: string;
  resolverFocusKey?: string;
}

export function createInitialState(): AppState {
  return {
    status: "Drop a Paperpile DOCX or choose a local file.",
  };
}

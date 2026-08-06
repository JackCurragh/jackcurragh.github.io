export type PdfRect = {
  x1: number;
  y1: number;
  x2: number;
  y2: number;
};

export type PaperpileReferenceId = string;

export interface PaperpileCitationAnnotation {
  pageIndex: number;
  rect: PdfRect;
  documentId: string;
  referenceIds: PaperpileReferenceId[];
  uri: string;
}

export interface PaperpileBibliographyAnnotation {
  pageIndex: number;
  rect: PdfRect;
  documentId: string;
  referenceId: PaperpileReferenceId;
  uri: string;
}

export interface BibliographyDestination {
  referenceId: PaperpileReferenceId;
  pageIndex: number;
  x?: number;
  y: number;
}

export interface InternalCitationLink {
  pageIndex: number;
  rect: PdfRect;
  referenceId: PaperpileReferenceId;
  destination: BibliographyDestination;
}

export interface ProcessingWarning {
  code: string;
  pageIndex?: number;
  message: string;
  severity: "info" | "warning" | "error";
}

export interface CitationProblem {
  code: string;
  problemKey: string;
  pageIndex: number;
  documentId: string;
  referenceIds: PaperpileReferenceId[];
  uri: string;
  message: string;
  evidenceText?: string;
  annotationRects: PdfRect[];
  previewDataUrl?: string;
  fragments: CitationProblemFragment[];
}

export interface CitationProblemFragment {
  index: number;
  text: string;
  rect: PdfRect;
}

export type AnnotationKind = "paperpile-citation" | "paperpile-bibliography" | "other-uri" | "non-uri";

export interface ClassifiedAnnotation {
  kind: AnnotationKind;
  pageIndex: number;
  rect: PdfRect;
  uri?: string;
  subtype?: string;
  documentId?: string;
  referenceIds?: PaperpileReferenceId[];
  referenceId?: PaperpileReferenceId;
}

export interface ProcessingReport {
  format: "docx";
  citationAnnotations: PaperpileCitationAnnotation[];
  bibliographyAnnotations: PaperpileBibliographyAnnotation[];
  bibliographyDestinations: BibliographyDestination[];
  internalLinks: InternalCitationLink[];
  warnings: ProcessingWarning[];
  problems: CitationProblem[];
  resolvedClusters: number;
  unresolvedClusters: number;
  figureTargets: string[];
  figureLinksCreated: number;
}

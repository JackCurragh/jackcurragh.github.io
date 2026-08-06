import type { ProcessingReport } from "../model/types";

export function summarizeProcessing(report: ProcessingReport): string {
  const citationGroupCount = new Set(
    report.citationAnnotations.map((citation) => `${citation.pageIndex}::${citation.documentId}::${citation.uri}`),
  ).size;
  const lines = [
    `DOCX Paperpile citation hyperlinks: ${report.citationAnnotations.length}`,
    `Paperpile citation groups: ${citationGroupCount}`,
    `DOCX bibliography hyperlinks: ${report.bibliographyAnnotations.length}`,
    `Bibliography bookmarks: ${report.bibliographyDestinations.length}`,
    `Resolved citation clusters: ${report.resolvedClusters}`,
    `Ambiguous citation clusters: ${report.unresolvedClusters}`,
    `Internal DOCX links created: ${report.internalLinks.length}`,
    `Figure and box targets: ${report.figureTargets.length}`,
    `Figure and box links created: ${report.figureLinksCreated}`,
  ];

  if (report.warnings.length > 0) lines.push(`Warnings: ${report.warnings.length}`);
  if (report.problems.length > 0) lines.push(`Resolver problems: ${report.problems.length}`);
  return lines.join("\n");
}

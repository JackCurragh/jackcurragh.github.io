export function collapseWhitespace(input: string): string {
  return input.replace(/\s+/g, " ").trim();
}

export function trimCitationPunctuation(input: string): string {
  return input.replace(/^[\s([{\u00a0]+|[\s)\]};,.:!?]+$/g, "").trim();
}

export function normalizeCitationText(input: string): string {
  return collapseWhitespace(input).replace(/\u00a0/g, " ");
}

export function splitSemicolonCluster(text: string): string[] {
  const segments = text.split(/\s*;\s*/g).map((segment) => trimCitationPunctuation(segment)).filter(Boolean);
  return segments;
}

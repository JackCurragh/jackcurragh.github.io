import type { PaperpileReferenceId } from "../model/types";

export type PaperpileUrlKind = "citation" | "bibliography";

export interface ParsedPaperpileUrl {
  kind: PaperpileUrlKind;
  documentId: string;
  referenceIds: PaperpileReferenceId[];
  uri: string;
}

const ALLOWED_HOSTS = new Set(["paperpile.com", "www.paperpile.com"]);

export function parsePaperpileUrl(uri: string): ParsedPaperpileUrl | undefined {
  let parsed: URL;
  try {
    parsed = new URL(uri);
  } catch {
    return undefined;
  }

  if (parsed.protocol !== "https:" || !ALLOWED_HOSTS.has(parsed.hostname)) {
    return undefined;
  }

  const pathSegments = parsed.pathname.split("/").filter(Boolean);
  if (pathSegments.length < 3) return undefined;

  const [scope, documentIdRaw, referenceSegmentRaw] = pathSegments;
  if (scope !== "c" && scope !== "b") return undefined;

  const documentId = safeDecode(pathSegments[1]);
  if (!documentId) return undefined;

  if (scope === "b") {
    if (pathSegments.length !== 3) return undefined;
    const referenceId = safeDecode(referenceSegmentRaw);
    if (!referenceId) return undefined;
    return {
      kind: "bibliography",
      documentId,
      referenceIds: [referenceId],
      uri,
    };
  }

  if (pathSegments.length !== 3) return undefined;
  const referenceIds = referenceSegmentRaw.split("+").map((part) => safeDecode(part)).filter((part): part is string => Boolean(part));
  if (referenceIds.length === 0 || referenceIds.some((ref) => ref.length === 0)) return undefined;
  return {
    kind: "citation",
    documentId,
    referenceIds,
    uri,
  };
}

function safeDecode(value: string): string | undefined {
  try {
    return decodeURIComponent(value);
  } catch {
    return undefined;
  }
}

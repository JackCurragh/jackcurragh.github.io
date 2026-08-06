import { strFromU8, strToU8, unzipSync, zipSync } from "fflate";
import type { CitationProblem, PdfRect, ProcessingReport } from "../model/types";
import { normalizeCitationText, splitSemicolonCluster } from "../utils/strings";
import { parsePaperpileUrl } from "../paperpile/urls";

const WORD_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main";
const REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";
const PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships";
const ZERO_RECT: PdfRect = { x1: 0, y1: 0, x2: 0, y2: 0 };

type Relationship = { id: string; target: string };
type ParsedLink = { element: Element; url: string; parsed: NonNullable<ReturnType<typeof parsePaperpileUrl>> };
type CitationLinkGroup = ParsedLink & { elements: Element[] };
type FigureTarget = { key: string; anchor: string; paragraph: Element };

export interface DocxProcessResult {
  report: ProcessingReport;
  outputBytes: Uint8Array;
}

export async function processDocx(file: File, progress: (message: string) => void): Promise<DocxProcessResult> {
  progress("Reading DOCX package locally…");
  const bytes = new Uint8Array(await file.arrayBuffer());
  const archive = unzipSync(bytes);
  const documentXml = requireEntry(archive, "word/document.xml");
  const relsXml = requireEntry(archive, "word/_rels/document.xml.rels");
  const document = parseXml(strFromU8(documentXml), "word/document.xml");
  const relationships = readRelationships(strFromU8(relsXml));
  const relationshipById = new Map(relationships.map((relationship) => [relationship.id, relationship.target]));
  const hyperlinks = [...document.getElementsByTagNameNS(WORD_NS, "hyperlink")];
  const citations = hyperlinks
    .map((element) => ({ element, url: relationshipById.get(element.getAttributeNS(REL_NS, "id") ?? "") }))
    .filter((item): item is { element: Element; url: string } => Boolean(item.url))
    .map((item) => ({ ...item, parsed: parsePaperpileUrl(item.url) }))
    .filter((item): item is { element: Element; url: string; parsed: NonNullable<ReturnType<typeof parsePaperpileUrl>> } => Boolean(item.parsed));

  const bibliography = citations.filter((item) => item.parsed.kind === "bibliography");
  const citationLinks = citations.filter((item) => item.parsed.kind === "citation");
  const citationGroups = groupAdjacentCitationLinks(citationLinks);
  const bibliographyIds = new Set(bibliography.map((item) => item.parsed.referenceIds[0]));
  const documentId = citations[0]?.parsed.documentId;
  const problems: CitationProblem[] = [];
  const internalLinks: Array<{ referenceId: string }> = [];
  let bookmarkId = nextBookmarkId(document);

  progress("Indexing figure and box captions…");
  const figureTargets = addFigureBookmarks(document, bookmarkId);
  bookmarkId += figureTargets.length;

  progress("Indexing bibliography bookmarks…");
  for (const item of bibliography) {
    const referenceId = item.parsed.referenceIds[0];
    const anchor = bookmarkName(referenceId);
    const start = document.createElementNS(WORD_NS, "w:bookmarkStart");
    start.setAttributeNS(WORD_NS, "w:id", String(bookmarkId));
    start.setAttributeNS(WORD_NS, "w:name", anchor);
    const end = document.createElementNS(WORD_NS, "w:bookmarkEnd");
    end.setAttributeNS(WORD_NS, "w:id", String(bookmarkId));
    bookmarkId += 1;
    unwrapHyperlink(item.element, [start, end]);
    internalLinks.push({ referenceId });
  }

  progress("Resolving Paperpile citation links…");
  for (const group of citationGroups) {
    const referenceIds = group.parsed.referenceIds;
    const missing = referenceIds.filter((referenceId) => !bibliographyIds.has(referenceId));
    const text = normalizeCitationText(group.elements.map((element) => element.textContent ?? "").join(""));
    const fragments = splitCitationText(text);
    const problemKey = `docx::${group.url}`;

    if (missing.length > 0) {
      problems.push(makeProblem(problemKey, documentId ?? group.parsed.documentId, referenceIds, group.url, "missing-bibliography", `No bibliography hyperlink exists for: ${missing.join(", ")}.`, text, fragments));
      for (const element of group.elements) unwrapHyperlink(element);
      continue;
    }

    if (referenceIds.length !== 1) {
      const rawSegments = splitRawCitationText(text, referenceIds.length);
      if (rawSegments.length === referenceIds.length) {
        splitCitationHyperlinks(group.elements, rawSegments, referenceIds);
        internalLinks.push(...referenceIds.map((referenceId) => ({ referenceId })));
      } else {
        problems.push(makeProblem(problemKey, documentId ?? group.parsed.documentId, referenceIds, group.url, "docx-multi-reference", "This DOCX hyperlink contains multiple Paperpile IDs, but its text did not split into the same number of citation segments.", text, fragments));
        for (const element of group.elements) unwrapHyperlink(element);
      }
      continue;
    }

    for (const element of group.elements) {
      element.removeAttributeNS(REL_NS, "id");
      element.setAttributeNS(WORD_NS, "w:anchor", bookmarkName(referenceIds[0]));
    }
    internalLinks.push({ referenceId: referenceIds[0] });
  }

  const figureLinkResult = linkFigureReferences(document, figureTargets);
  const figureWarnings = figureLinkResult.warnings.map((warning) => ({
    code: "missing-figure-target",
    severity: "warning" as const,
    message: warning.message,
  }));

  const outputXml = new XMLSerializer().serializeToString(document);
  const outputArchive = { ...archive, "word/document.xml": strToU8(outputXml) };
  const report: ProcessingReport = {
    format: "docx",
    citationAnnotations: citationLinks.map((item) => ({ pageIndex: -1, rect: ZERO_RECT, documentId: item.parsed.documentId, referenceIds: item.parsed.referenceIds, uri: item.url })),
    bibliographyAnnotations: bibliography.map((item) => ({ pageIndex: -1, rect: ZERO_RECT, documentId: item.parsed.documentId, referenceId: item.parsed.referenceIds[0], uri: item.url })),
    bibliographyDestinations: bibliography.map((item) => ({ referenceId: item.parsed.referenceIds[0], pageIndex: -1, y: 0 })),
    internalLinks: internalLinks.map((item) => ({ pageIndex: -1, rect: ZERO_RECT, referenceId: item.referenceId, destination: { referenceId: item.referenceId, pageIndex: -1, y: 0 } })),
    warnings: [
      ...problems.map((problem) => ({ code: problem.code, severity: "warning" as const, message: problem.message })),
      ...figureWarnings,
    ],
    problems,
    resolvedClusters: citationGroups.length - problems.length,
    unresolvedClusters: problems.length,
    figureTargets: figureTargets.map((target) => target.key),
    figureLinksCreated: figureLinkResult.linksCreated,
  };

  progress("Building rewritten DOCX…");
  return { report, outputBytes: zipSync(outputArchive) };
}

function nextBookmarkId(document: XMLDocument): number {
  let maximum = 0;
  for (const bookmark of [...document.getElementsByTagNameNS(WORD_NS, "bookmarkStart")]) {
    const id = Number.parseInt(bookmark.getAttributeNS(WORD_NS, "id") ?? "", 10);
    if (Number.isFinite(id)) maximum = Math.max(maximum, id);
  }
  return maximum + 1;
}

function readRelationships(xml: string): Relationship[] {
  const document = parseXml(xml, "word/_rels/document.xml.rels");
  return [...document.getElementsByTagNameNS(PACKAGE_REL_NS, "Relationship")].map((element) => ({
    id: element.getAttribute("Id") ?? "",
    target: element.getAttribute("Target") ?? "",
  })).filter((relationship) => relationship.id && relationship.target);
}

function addFigureBookmarks(document: XMLDocument, firstBookmarkId: number): FigureTarget[] {
  const targets: FigureTarget[] = [];
  let bookmarkId = firstBookmarkId;
  for (const paragraph of [...document.getElementsByTagNameNS(WORD_NS, "p")]) {
    const key = captionKey(paragraph.textContent ?? "");
    if (!key || targets.some((target) => target.key === key)) continue;
    const anchor = key;
    const start = document.createElementNS(WORD_NS, "w:bookmarkStart");
    start.setAttributeNS(WORD_NS, "w:id", String(bookmarkId));
    start.setAttributeNS(WORD_NS, "w:name", anchor);
    const end = document.createElementNS(WORD_NS, "w:bookmarkEnd");
    end.setAttributeNS(WORD_NS, "w:id", String(bookmarkId));
    bookmarkId += 1;

    const firstContent = [...paragraph.childNodes].find((node) => node.nodeType === 1 && (node as Element).localName !== "pPr");
    paragraph.insertBefore(start, firstContent ?? null);
    paragraph.append(end);
    targets.push({ key, anchor, paragraph });
  }
  return targets;
}

function linkFigureReferences(document: XMLDocument, targets: FigureTarget[]): { linksCreated: number; warnings: Array<{ key: string; text: string; message: string }> } {
  const targetByKey = new Map(targets.map((target) => [target.key, target]));
  const captionParagraphs = new Set(targets.map((target) => target.paragraph));
  const warnings: Array<{ key: string; text: string; message: string }> = [];
  let linksCreated = 0;

  for (const paragraph of [...document.getElementsByTagNameNS(WORD_NS, "p")]) {
    if (captionParagraphs.has(paragraph)) continue;
    for (const run of [...paragraph.getElementsByTagNameNS(WORD_NS, "r")]) {
      if (hasHyperlinkAncestor(run)) continue;
      const textNodes = [...run.getElementsByTagNameNS(WORD_NS, "t")];
      if (textNodes.length !== 1) continue;
      const text = textNodes[0].textContent ?? "";
      const matches = [...text.matchAll(/\b(Fig(?:ure)?\.?|Box)\s+(\d+)([a-z])?\b/gi)];
      if (matches.length === 0) continue;

      const parent = run.parentNode;
      if (!parent) continue;
      let cursor = 0;
      for (const match of matches) {
        const before = text.slice(cursor, match.index);
        if (before) parent.insertBefore(cloneRunWithText(run, before), run);
        const key = `${match[1].toLowerCase().startsWith("box") ? "box" : "figure"}-${match[2]}`;
        const target = targetByKey.get(key);
        const matchedText = match[0];
        if (target) {
          const hyperlink = document.createElementNS(WORD_NS, "w:hyperlink");
          hyperlink.setAttributeNS(WORD_NS, "w:anchor", target.anchor);
          hyperlink.setAttributeNS(WORD_NS, "w:history", "1");
          hyperlink.append(cloneRunWithText(run, matchedText));
          parent.insertBefore(hyperlink, run);
          linksCreated += 1;
        } else {
          parent.insertBefore(cloneRunWithText(run, matchedText), run);
          warnings.push({ key: `docx::${key}`, text: matchedText, message: `No ${key.startsWith("box") ? "box" : "figure"} caption was found for ${matchedText}.` });
        }
        cursor = (match.index ?? 0) + matchedText.length;
      }
      const after = text.slice(cursor);
      if (after) parent.insertBefore(cloneRunWithText(run, after), run);
      parent.removeChild(run);
    }
  }
  return { linksCreated, warnings };
}

function captionKey(text: string): string | undefined {
  const match = text.match(/^\s*(Figure|Fig\.?|Box)\s+(\d+)(?:[a-z])?\b/i);
  if (!match) return undefined;
  return `${match[1].toLowerCase().startsWith("box") ? "box" : "figure"}-${match[2]}`;
}

function hasHyperlinkAncestor(element: Element): boolean {
  let parent = element.parentElement;
  while (parent) {
    if (parent.localName === "hyperlink") return true;
    parent = parent.parentElement;
  }
  return false;
}

function cloneRunWithText(run: Element, text: string): Element {
  const clone = run.ownerDocument!.createElementNS(WORD_NS, "w:r");
  const properties = run.getElementsByTagNameNS(WORD_NS, "rPr")[0];
  if (properties) clone.append(properties.cloneNode(true));
  const textNode = run.ownerDocument!.createElementNS(WORD_NS, "w:t");
  textNode.textContent = text;
  if (/^\s|\s$/.test(text)) textNode.setAttributeNS("http://www.w3.org/XML/1998/namespace", "xml:space", "preserve");
  clone.append(textNode);
  return clone;
}

function groupAdjacentCitationLinks(links: ParsedLink[]): CitationLinkGroup[] {
  const groups: CitationLinkGroup[] = [];
  for (const link of links) {
    const previous = groups[groups.length - 1];
    if (previous && previous.url === link.url && previous.element.parentNode === link.element.parentNode && areAdjacent(previous.elements[previous.elements.length - 1], link.element)) {
      previous.elements.push(link.element);
      continue;
    }
    groups.push({ ...link, elements: [link.element] });
  }
  return groups;
}

function areAdjacent(first: Element, second: Element): boolean {
  let node = first.nextSibling;
  while (node && node !== second) {
    if (node.textContent?.trim()) return false;
    node = node.nextSibling;
  }
  return node === second;
}

function unwrapHyperlink(element: Element, wrappers: Element[] = []): void {
  const parent = element.parentNode;
  if (!parent) return;
  const fragment = element.ownerDocument!.createDocumentFragment();
  for (const wrapper of wrappers) fragment.append(wrapper);
  while (element.firstChild) fragment.append(element.firstChild);
  for (const child of [...fragment.childNodes]) parent.insertBefore(child, element);
  parent.removeChild(element);
}

function splitCitationText(text: string): Array<{ index: number; text: string; rect: PdfRect }> {
  return splitSemicolonCluster(text).map((segment, index) => ({ index, text: segment, rect: ZERO_RECT }));
}

function splitRawCitationText(text: string, expectedCount: number): string[] {
  const segments: string[] = [];
  let start = 0;
  for (let index = 0; index < text.length; index += 1) {
    if (text[index] !== ";") continue;
    segments.push(text.slice(start, index + 1));
    start = index + 1;
  }
  if (start < text.length) segments.push(text.slice(start));
  return segments.length === expectedCount ? segments : [];
}

function splitCitationHyperlinks(elements: Element[], segments: string[], referenceIds: string[]): void {
  const element = elements[0];
  const parent = element.parentNode;
  const document = element.ownerDocument;
  if (!parent || !document) return;
  const runs = [...element.getElementsByTagNameNS(WORD_NS, "r")];
  const template = runs[0]?.cloneNode(true) as Element | undefined;
  if (!template) return;

  for (let index = 0; index < segments.length; index += 1) {
    const hyperlink = element.cloneNode(false) as Element;
    hyperlink.removeAttributeNS(REL_NS, "id");
    hyperlink.setAttributeNS(WORD_NS, "w:anchor", bookmarkName(referenceIds[index]));
    const run = template.cloneNode(true) as Element;
    const textNodes = [...run.getElementsByTagNameNS(WORD_NS, "t")];
    if (textNodes.length === 0) {
      const textNode = document.createElementNS(WORD_NS, "w:t");
      run.append(textNode);
      textNodes.push(textNode);
    }
    textNodes[0].textContent = segments[index];
    for (const extraTextNode of textNodes.slice(1)) extraTextNode.remove();
    if (/^\s|\s$/.test(segments[index])) textNodes[0].setAttributeNS("http://www.w3.org/XML/1998/namespace", "xml:space", "preserve");
    hyperlink.append(run);
    parent.insertBefore(hyperlink, element);
  }
  for (const original of elements) parent.removeChild(original);
}

function makeProblem(problemKey: string, documentId: string, referenceIds: string[], uri: string, code: string, message: string, evidenceText: string, fragments: Array<{ index: number; text: string; rect: PdfRect }>): CitationProblem {
  return { code, problemKey, pageIndex: -1, documentId, referenceIds, uri, message, evidenceText, annotationRects: [], fragments };
}

function bookmarkName(referenceId: string): string {
  return `citation-${referenceId.replace(/[^A-Za-z0-9_.-]/g, "-")}`;
}

function requireEntry(archive: Record<string, Uint8Array>, path: string): Uint8Array {
  const entry = archive[path];
  if (!entry) throw new Error(`The DOCX package is missing ${path}.`);
  return entry;
}

function parseXml(xml: string, source: string): XMLDocument {
  const document = new DOMParser().parseFromString(xml, "application/xml");
  if (document.querySelector("parsererror")) throw new Error(`Could not parse ${source}.`);
  return document;
}

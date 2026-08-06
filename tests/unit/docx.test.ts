import { strFromU8, strToU8, unzipSync, zipSync } from "fflate";
import { describe, expect, it } from "vitest";
import { processDocx } from "../../src/docx/process";

const documentXml = `<?xml version="1.0" encoding="UTF-8"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <w:body>
    <w:p><w:hyperlink r:id="rId1"><w:r><w:t>(Smith 2020)</w:t></w:r></w:hyperlink></w:p>
    <w:p><w:hyperlink r:id="rId2"><w:r><w:t>Smith, J. (2020). A paper.</w:t></w:r></w:hyperlink></w:p>
  </w:body>
</w:document>`;

const relationshipsXml = `<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" Target="https://paperpile.com/c/doc1/ref1" TargetMode="External"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink" Target="https://paperpile.com/b/doc1/ref1" TargetMode="External"/>
</Relationships>`;

describe("processDocx", () => {
  it("converts Paperpile citation and bibliography hyperlinks to internal DOCX anchors", async () => {
    const bytes = zipSync({
      "word/document.xml": strToU8(documentXml),
      "word/_rels/document.xml.rels": strToU8(relationshipsXml),
    });
    const input = { arrayBuffer: async () => bytes.buffer } as unknown as File;

    const result = await processDocx(input, () => undefined);
    const output = strFromU8(unzipSync(result.outputBytes)["word/document.xml"]);

    expect(result.report.problems).toHaveLength(0);
    expect(result.report.resolvedClusters).toBe(1);
    expect(result.report.internalLinks).toHaveLength(2);
    expect(output).toContain("w:anchor=\"citation-ref1\"");
    expect(output).toContain("w:bookmarkStart");
    expect(output).toContain("w:name=\"citation-ref1\"");
    expect(output).not.toContain("r:id=\"rId1\"");
    expect(output).not.toContain("r:id=\"rId2\"");
  });

  it("splits semicolon-separated multi-reference hyperlinks safely", async () => {
    const multiDocument = documentXml
      .replace("(Smith 2020)", "(Smith 2020; Jones 2021)")
      .replace("rId1", "rId3")
      .replace("</w:body>", "<w:p><w:hyperlink r:id=\"rId4\"><w:r><w:t>Second paper.</w:t></w:r></w:hyperlink></w:p></w:body>");
    const multiRelationships = relationshipsXml
      .replace("Id=\"rId1\"", "Id=\"rId3\"")
      .replace("/ref1\" TargetMode", "/ref1+ref2\" TargetMode")
      .replace("</Relationships>", "<Relationship Id=\"rId4\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink\" Target=\"https://paperpile.com/b/doc1/ref2\" TargetMode=\"External\"/></Relationships>");
    const bytes = zipSync({
      "word/document.xml": strToU8(multiDocument),
      "word/_rels/document.xml.rels": strToU8(multiRelationships),
    });
    const result = await processDocx({ arrayBuffer: async () => bytes.buffer } as unknown as File, () => undefined);

    const output = strFromU8(unzipSync(result.outputBytes)["word/document.xml"]);
    expect(result.report.problems).toHaveLength(0);
    expect(result.report.resolvedClusters).toBe(1);
    expect(result.report.internalLinks).toHaveLength(4);
    expect(output.match(/w:anchor="citation-ref[12]"/g)).toHaveLength(2);
  });

  it("links figure and box references to their captions", async () => {
    const figureDocument = `<?xml version="1.0" encoding="UTF-8"?>
      <w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
        <w:body>
          <w:p><w:r><w:t>Figure 1. Experimental setup</w:t></w:r></w:p>
          <w:p><w:r><w:t>See Fig. 1 and Box 2 for details.</w:t></w:r></w:p>
          <w:p><w:r><w:t>Box 2. Analysis notes</w:t></w:r></w:p>
        </w:body>
      </w:document>`;
    const emptyRelationships = `<?xml version="1.0" encoding="UTF-8"?>
      <Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>`;
    const bytes = zipSync({
      "word/document.xml": strToU8(figureDocument),
      "word/_rels/document.xml.rels": strToU8(emptyRelationships),
    });

    const result = await processDocx({ arrayBuffer: async () => bytes.buffer } as unknown as File, () => undefined);
    const output = strFromU8(unzipSync(result.outputBytes)["word/document.xml"]);

    expect(result.report.problems).toHaveLength(0);
    expect(result.report.figureTargets).toEqual(["figure-1", "box-2"]);
    expect(result.report.figureLinksCreated).toBe(2);
    expect(output).toContain('w:name="figure-1"');
    expect(output).toContain('w:name="box-2"');
    expect(output).toContain('w:anchor="figure-1"');
    expect(output).toContain('w:anchor="box-2"');
  });
});

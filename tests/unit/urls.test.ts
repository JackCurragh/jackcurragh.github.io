import { describe, expect, it } from "vitest";
import { parsePaperpileUrl } from "../../src/paperpile/urls";

describe("parsePaperpileUrl", () => {
  it("parses a single bibliography URL", () => {
    expect(parsePaperpileUrl("https://paperpile.com/b/Eumm7u9g/wnCR")).toEqual({
      kind: "bibliography",
      documentId: "Eumm7u9g",
      referenceIds: ["wnCR"],
      uri: "https://paperpile.com/b/Eumm7u9g/wnCR",
    });
  });

  it("parses a citation cluster", () => {
    expect(parsePaperpileUrl("https://www.paperpile.com/c/Eumm7u9g/wnCR+0eWX+E6nr")).toEqual({
      kind: "citation",
      documentId: "Eumm7u9g",
      referenceIds: ["wnCR", "0eWX", "E6nr"],
      uri: "https://www.paperpile.com/c/Eumm7u9g/wnCR+0eWX+E6nr",
    });
  });

  it("rejects non-Paperpile URLs", () => {
    expect(parsePaperpileUrl("https://example.com/c/Eumm7u9g/wnCR")).toBeUndefined();
  });
});

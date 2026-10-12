# Java Straw

Java Straw is an in-process Java reader for indexed access to Juicebox `.hic`
files. It supports local files and URLs, V6–V9 and consolidated V10 files,
sparse block/iterator access, BP and FRAG resolutions, normalization vectors,
expected values, Pearson matrices, eigenvectors, and dynamic resolutions.

## Usage

```java
Dataset dataset = HiCFileTools.extractDatasetForCLT("sample.hic", false, true, false);
Chromosome chr2 = dataset.getChromosomeHandler().getChromosomeFromName("2");
Chromosome chr1 = dataset.getChromosomeHandler().getChromosomeFromName("1");

// X is chr2 and Y is chr1 because that is the order requested here.
Matrix matrix = dataset.getMatrix(chr2, chr1);
MatrixZoomData zoom = matrix.getZoomData(new HiCZoom(10_000));
List<Block> blocks = zoom.getNormalizedBlocksOverlapping(
    0, 0, 100, 100, NormalizationHandler.NONE, false);

for (Block block : blocks) {
    for (ContactRecord record : block.getContactRecords()) {
        System.out.println(record.getBinX() + "\t" + record.getBinY() + "\t" +
                           record.getCounts());
    }
}
```

See [`src/javastraw/AnnotatedExample.java`](src/javastraw/AnnotatedExample.java)
for iterator and window-query examples.

## Axis and count contracts

For inter-chromosomal matrices, `Dataset.getMatrix(first, second)` preserves the
caller’s order: X coordinates and `MatrixZoomData.getChr1()` refer to `first`,
while Y and `getChr2()` refer to `second`. Canonical on-disk ordering remains an
internal detail. This applies to block queries, direct/normalized iterators, and
dynamic-resolution views.

`ContactRecord.getCounts()` remains `float` for source compatibility. For a V10
integer-valued record, use `ContactRecord.getExactCount()`, which returns
`BigInteger` and preserves the complete unsigned range `0..2^64-1`, including
derived-resolution aggregation. V10 score records remain floating point.

## V10 implementation

The V10 reader implements the final block-index format, not the withdrawn
page-based draft. Each non-empty logical block is located by its exact `H10I`
version-2 entry and decompressed from one independent `H10B` Zstandard frame.
Regional reads compute exact block-number ranges (including the nearest rotated
cis distance band), binary-search the cached index, and decode positions and
values in one streaming pass. Parsed matrix descriptors, block indexes,
decompressed blocks, and per-thread Zstandard decompressors are reused across
queries. Physically adjacent requested blocks may share one range read, but an
unrequested block is never decompressed.

The standard v10 resolution pyramid is a writer default, not a reader constraint.
Java Straw follows each resolution's declared materialized/derived mode and direct
materialized source, including custom base resolutions and derivation patterns.

For V10 trans O/E, `MatrixZoomData.getAverageCount()` supplies the raw full-matrix
mean `(sumCountsOrScores / nBins1) / nBins2`, including for score matrices.
The bin counts follow V10 geometry (partial BP bins and FRAG bins) at the
requested logical resolution, including derived resolutions. Existing trans
O/E extraction applies normalization to the numerator and uses this raw mean
as its denominator, as with V9; no cis expected vectors are needed.
`getExpectedValueForTrans()` marks zero or non-finite V10 means unavailable;
finite negative score means remain valid for arithmetic O/E. Older formats
retain Java Straw's existing fallback behavior.

## Build

The repository includes its compile-time JAR dependencies under `lib/`:

```sh
mkdir -p build/classes
find src -name '*.java' -print0 | \
  xargs -0 javac -cp 'lib/broadinstitute/*:lib/general/*' -d build/classes
jar --create --file build/java-straw.jar -C build/classes .
```

On Windows, replace the classpath `:` separator with `;`. IntelliJ users can
also use the checked-in project metadata and artifact configuration.

Slice and HBS creation are intentionally C++-only and are not implemented in
Java.

## Support

Use the [Aiden Lab forum](https://aidenlab.org/forum.html) for questions and the
repository issue tracker for reproducible bugs or feature requests.


## Stored scaled All in V10

Files produced by current hictools V10 writers may append `All` at any header
index. Its useful overview is stored at the smallest advertised BP resolution
in synthetic coordinates. No All-to-real matrices are required; absent pairs
query as zero. Real chromosome queries and the global resolution list retain
their ordinary semantics. Query raw `NONE` data for All; it has no normalization
vectors. Coarser All levels are valid compact aggregates.

Version-one metadata uses existing attributes: `hicAll.version`,
`hicAll.finestResolution`, `hicAll.sourceResolution`, `hicAll.groupFactor`,
`hicAll.binCount`, and `hicAll.chromosomeBins`. The last value lists padded
source-bin counts of real chromosomes in header order. Physical overview bin
size is source resolution times group factor; raw queries use the finest
synthetic resolution. This convention changes no binary fields or file version.

`DatasetReaderV10.getHeader().allOverview` exposes the validated layout.
All filtering uses its name at any table index; exact v10 chromosome names are preserved.

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

#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes=${TMPDIR:-/tmp}/java-straw-test-$$
trap 'rm -rf "$classes"' EXIT HUP INT TERM
mkdir -p "$classes"
find "$root/src" "$root/tests" -name '*.java' -print0 | \
  xargs -0 javac -proc:none -cp "$root/lib/broadinstitute/*:$root/lib/general/*" -d "$classes"
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10UnsignedRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10BlockIndexRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10BlockDecoderRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10ResolutionPolicyRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" CallerOrientationRegression "$root/data/inter.hic"
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" TransExpectedRegression \
  "$root/tests/data/trans-v10.hic" "$root/tests/data/trans-score-v10.hic"

java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" NormalizationValidityRegression "$root/tests/data/cis-v10.hic"
if [ "${JAVA_STRAW_TEST_HTTP:-0}" = 1 ]; then
  java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" HttpRangeRegression "$root/tests/data/cis-v10.hic"
fi

java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" ScaledAllRegression "$root/tests/data/pseudo-all.hic" "$root/tests/data/kilobase-all.hic"
java -ea -Xmx128m -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10ReviewRegression "$root/tests/data"

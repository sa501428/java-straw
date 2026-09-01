#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
classes=${TMPDIR:-/tmp}/java-straw-test-$$
trap 'rm -rf "$classes"' EXIT HUP INT TERM
mkdir -p "$classes"
find "$root/src" "$root/tests" -name '*.java' -print0 | \
  xargs -0 javac -cp "$root/lib/broadinstitute/*:$root/lib/general/*" -d "$classes"
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10UnsignedRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10BlockIndexRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" V10BlockDecoderRegression
java -ea -cp "$classes:$root/lib/broadinstitute/*:$root/lib/general/*" CallerOrientationRegression "$root/data/inter.hic"

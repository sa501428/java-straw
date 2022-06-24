--------------
About Java Straw
--------------
Java Straw is a library for quickly reading data from the HiC file format, to be used directly in code without dumping
data to a local file.

Here's a detailed example of how to use Java-Straw:
[AnnotatedExample.java](https://github.com/sa501428/java-straw/blob/master/src/javastraw/AnnotatedExample.java)

Ongoing development work is carried out by <a href="https://github.com/sa501428">Muhammad Saad Shamim</a>. Past
contributors include <a href="https://github.com/nchernia">Neva C. Durand</a>, <a href="https://github.com/suhas-rao">
Suhas Rao</a>, and <a href="https://github.com/jrobinso">Jim Robinson</a>.

--------------
Questions?
--------------

For FAQs, or for asking new questions, please see our forum: <a href="https://aidenlab.org/forum.html">
aidenlab.org/forum.html</a>.

--------------
Compiling Jar from IntelliJ
--------------

Precompiled jars are available under [Github Releases](https://github.com/sa501428/java-straw/releases). We recommend
using the latest version.

You can build jars from source code using IntelliJ IDEA (Community edition - free)

To set up in IDEA, have the Java SDK installed then you'll point to it (IntelliJ has lots of documentation on this sort
of thing).

* Then go to `VCS` -> `checkout from version control`.
* You'll need to be sure `*.sizes` is included as a file to be copied over to the class files. Set this up via
  IntelliJ `Preferences` -> `Compiler`. Add `?*.sizes` to the list of `Resource Patterns`.
* One last note: be sure to `Commit and Push` when you commit files, it's hidden in the dropdown menu button in the
  commit window.
* Compiling the code should automatically build the jar under out/artifacts/java_straw_jar

-------------
Documentation
-------------
We have documentation for how to use Java Straw at
https://github.com/aidenlab/java-straw/wiki.


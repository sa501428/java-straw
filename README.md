--------------
About Java Straw
--------------
Java Straw is a library for quickly reading data from the HiC file format, to be used directly in code without dumping
data to a local file.

Ongoing development work is carried out by <a href="https://github.com/sa501428">Muhammad Saad Shamim</a>. Past
contributors include <a href="https://github.com/nchernia">Neva C. Durand</a>, <a href="https://github.com/suhas-rao">
Suhas Rao</a>, and <a href="https://github.com/jrobinso">Jim Robinson</a>.

--------------
Questions?
--------------

For FAQs, or for asking new questions, please see our forum: <a href="https://aidenlab.org/forum.html">aidenlab.org/forum.html</a>.

--------------
IntelliJ Setup
--------------

Use IntelliJ IDEA (Community edition - free)

To set up in IDEA, have the Java SDK installed
then you'll point to it (IntelliJ has lots of documentation on this sort of thing).

* Then go to `VCS` -> `checkout from version control`.
* You'll need to do is be sure `*.sizes` is included as a file to be copied over to the class files.
Set this up via IntelliJ `Preferences` -> `Compiler`. Add `?*.sizes` to the list of `Resource Patterns`.
* While there, also go to `Java Compiler` and put this into additional command line options: `-Xlint:all -target 1.7`
The former turns on all warnings, the latter gives some flexibility since some people haven't updated Java to 1.8 yet.
* Then go to `Run` -> `Edit Configurations`.
* With the `+` sign, add `Application`.
* You'll create two of these, one for the GUI (call it Juicebox GUI or whatever you want, really) and one for the CLT.
* Set the main class by clicking the little `...` button next to the text box for main class

        HiCTools.java is the main method class for the analysis/CLT portion.

* For the CLT use

        -Xmx2000m

* Note that the `Xmx2000m` flag sets the maximum memory heap size to 2GB.
Depending on your computer you might want more or less.
Some tools will break if there's not enough memory and the file is too large,
but don't worry about that for development; 2GB should be fine.
* One last note: be sure to `Commit and Push` when you commit files, it's hidden in the dropdown menu button in the
commit window.

-------------
Documentation
-------------
We have documentation for how to use Java Straw at
https://github.com/aidenlab/java-straw/wiki.

--------------------------------
Compiling Jars from Source Files
--------------------------------
1. You should have Java 1.8 JDK and Apache Ant installed on your system. See
   below for more information.
2. Go to the folder containing the Juicebox source files and edit the
   javastraw.properties file with the proper Java JDK Address.
3. Open the command line, navigate to the folder containing the build.xml file
   and type
     ant
   The process should take no more than a minute to build on most machines.
4. The jars are written to the directory out/.  You can change this by editing
   the build.xml file.

* Installing Java 1.8 JDK

For Windows/Mac/Linux, the Java 1.8 JDK can be installed from here:
https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
(Alternative) For Ubuntu/LinuxMint
https://tecadmin.net/install-oracle-java-8-jdk-8-ubuntu-via-ppa/

* Installing Apache Ant
Mac
  Ant should be installed on most Macs. To verify installation via the command
  prompt, type
    ant -version
  If Ant is not on your Mac, install it via homebrew. At the command prompt, type
    brew update
    brew install ant
  You may need to install Homebrew (https://brew.sh/) on your machine
  See the following Stackoverflow post for more details:
  https://stackoverflow.com/questions/3222804/how-can-i-install-apache-ant-on-mac-os-x

Windows
  Installing Ant requires some minor changes to your system environment. Follow
  the instructions in this article:
  https://www.nczonline.net/blog/2012/04/12/how-to-install-apache-ant-on-windows/

Linux
  In the command prompt, type
    sudo apt-get install ant
  or
    sudo yum install ant
  depending on your package installer

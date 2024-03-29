/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2011-2020 Broad Institute, Aiden Lab, Rice University, Baylor College of Medicine
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 *  THE SOFTWARE.
 */

package javastraw.tools;

import javastraw.StrawGlobals;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.List;

/**
 * Created by muhammadsaadshamim on 9/29/15.
 */
public class UNIXTools {
    public static String extractElement(String str, int i) {
        String[] strSplit = str.split("\t");
        return strSplit[strSplit.length - i];
    }

    public static String makeDir(String path) {
        File file = new File(path);
        if (!file.isDirectory()) {
            file.mkdir();
        }
        return path;
    }

    public static File makeDir(File folder) {
        if (!folder.isDirectory()) {
            folder.mkdir();
        }
        return folder;
    }

    public static void redirectOutput(List<String> command, String outputFilePath, boolean printVerbose) {
        String output = executeComplexCommand(command, printVerbose);
        File outputFile = new File(outputFilePath);
        try {
            outputFile.createNewFile();
            PrintWriter writer = new PrintWriter(outputFile);
            writer.print(output);
            writer.close();
        } catch (Exception e) {
            System.err.println("Unable to write command line output to file: " + outputFilePath);
            e.printStackTrace();
        }
    }

    public static String executeSimpleCommand(String command) {
        StringBuilder output = new StringBuilder();
        try {
            Process p = Runtime.getRuntime().exec(command);
            p.waitFor();
            BufferedReader reader =
                    new BufferedReader(new InputStreamReader(p.getInputStream()), StrawGlobals.bufferSize);

            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return output.toString();
    }

    private static String executeComplexCommand(List<String> command, boolean printVerboseComments) {
        StringBuilder output = new StringBuilder();

        //System.out.println(System.getenv());

        ProcessBuilder b = new ProcessBuilder(command);
        //Map<String, String> env = b.environment();
        //System.out.println(env);

        Process p;
        try {
            //p = Runtime.getRuntime().exec(command);
            p = b.redirectErrorStream(true).start();

            if (printVerboseComments) {
                System.out.println("Command exec " + p.waitFor());
            } else {
                p.waitFor();
            }

            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()), StrawGlobals.bufferSize);

            String line;
            while ((line = reader.readLine()) != null) {
                output.append(line).append("\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return output.toString();
    }
}

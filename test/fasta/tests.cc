#include <gtest/gtest.h>
#include <iostream>
#include "fasta/fastareader.h"
#include <string>
#include <bits/stdc++.h>


std::string referencepath = "../../../test/data/reference.fa";
std::string referenceindexpath = "../../../test/data/reference.fa.fai";

TEST(FilePathTest, getTestFasta1)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    fastaReader.exitIfNoFilePath();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",0,10);
    ASSERT_EQ("CGAGCCGAAC", seqrefString);
}

TEST(FilePathTest, getTestFasta2)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr2",0,10);
    ASSERT_EQ("GTCCTCGGTG", seqrefString);
}

TEST(FilePathTest, getTestFasta3)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",50,50+10);
    ASSERT_EQ("TAAAGGGGTC", seqrefString);
}

TEST(FilePathTest, getTestFasta4)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr2",75,75+10);
    ASSERT_EQ("GCCCCAGCTC", seqrefString);
}

TEST(FilePathTest, getTestFasta5)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",0,0);
    ASSERT_EQ("", seqrefString);
}

TEST(FilePathTest, getTestFasta6)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",50,50);
    ASSERT_EQ("", seqrefString);
}

TEST(FilePathTest, getTestFasta7)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",0,1);
    ASSERT_EQ("C", seqrefString);
}

TEST(FilePathTest, getTestFasta8)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",50,51);
    ASSERT_EQ("T", seqrefString);
}

TEST(FilePathTest, getTestFasta9)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr1",51,52);
    ASSERT_EQ("A", seqrefString);
}

TEST(FilePathTest, getTestFasta10)
{
    FastaReader fastaReader;
    fastaReader.setFilePath(referencepath);
    fastaReader.setIndexFilePath(referenceindexpath);
    fastaReader.initialize();
    std::string seqrefString = fastaReader.getSeqbyPosition("chr2",75,76);
    ASSERT_EQ("G", seqrefString);
}
# Bare Minimum Implementatioin of Nextflow DSL2 on Sumner

## Objective

Nextflow is a wrapper for working with linux containers. DSL2 stands for
domain specific language number 2. This is because the syntax from DSL1 has changed.
This means we are standardizing to DSL2 and will update old scripts to
meet its requirements.
<hr>

## Key Take Aways

This toy example will help you get up and started on the key structure that the
NF-Core uses. Such as how to configure your parameters, start a process, get an output,
push the output to the next process. There are much more advanced use cases
than this example. Once you understand the basics here you can move on to see our RNASeq Pipeline that utilizes some advanced features and structures.
<hr>

## Primary Resources

<ul>
<li><a href="https://www.nextflow.io/docs/latest/dsl2.html">Link to DSL2 Documentation</a></li>
<li><a href="https://github.com/nf-core/hic">Link to HiC NF-Core</a></li>
<li><a href="https://github.com/nf-core/rnafusion">Link to RNA Fusion NF-Core</a></li>
</ul
<hr>

## What This Script Does

In short, this script runs two processes it will trim raw Paired End RNASeqs that
are found in the $$$ directory on Sumner. It will then push the trimmed sequences to
RSEM process to quantify the trimmed reads. The expected output can be found in $$$ folder.
<hr>

## Other Things to Note

decision making logic (some do specify some do not) need to figure out dsl2 RNA fusion and hic pipeline.
srun -p compute -t 02:00:00 --mem 2G --pty bash

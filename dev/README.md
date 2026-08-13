# dev/ — author tooling, not part of the pipeline

Nothing in here is run by Snakemake, referenced by a rule, or needed to use the
pipeline. These are ad-hoc analysis scripts kept because they answer questions
the pipeline itself cannot.

## `diag_capture.sh`

Attributes every hap.py false negative to either an **uncovered** or a
**covered** region, by deriving the actually-callable region from read depth in
the recalibrated BAM.

This distinguishes the two very different explanations for a mediocre recall
number:

- **Capture-driven** — the truth variants sit in target regions the capture
  never covered. The caller is at its ceiling; neither more depth nor a
  different BED will recover them.
- **Pipeline-driven** — the regions were covered and the pipeline still missed
  them, so alignment, calling or filtering is at fault.

Run it after a full run plus the `benchmark` target; see the usage block at the
top of the script. It needs `samtools` and `bedtools` from the baked conda envs,
so it must run inside the project image.

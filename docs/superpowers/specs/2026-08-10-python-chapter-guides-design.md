# Python Chapter Guides Design

## Goal

Make every existing chapter under `python/` understandable to a reader who
does not have the accompanying book. The documentation should provide enough
scientific and mathematical context to understand what each program models,
why the experiment is useful, and how to interpret its output.

The guides are companions to the implementations, not a replacement textbook.
They will use original explanations of the underlying concepts and will not
reproduce the book's prose, exercises, extended derivations, or illustrations.
The supplied PDF may be used to verify topic alignment and terminology.

## Documentation Location

Add a `README.md` to each of the 38 existing chapter directories under
`python/`. Keeping the explanation next to the implementations makes it visible
when browsing on GitHub and keeps links to example directories short and
stable.

Add `python/README.md` as the entry point. It will contain:

- a short explanation of how the Python material is organized;
- a linked index of all documented chapter directories;
- a suggested learning path grouped into single-neuron models, single-neuron
  dynamics, neuronal communication, network synchronization and rhythms, and
  plasticity;
- a clear note that book Chapters 2 and 6 do not have Python example
  directories in this repository.

Add one link from the root `README.md` to this index. Do not create empty
directories for chapters without code.

## Chapter Numbering and Titles

The filesystem is the source of truth for documentation placement. Each guide
will identify both its repository directory and the corresponding scientific
topic. Where a directory number or title differs from the book's table of
contents, the guide will state the relationship explicitly rather than imply an
exact chapter match.

This is especially important for the first directory, whose name describes the
first modeled neuron although the book's numbered Chapter 1 is introductory
vocabulary, and for directories whose historical names differ from the final
book title.

## Standard Chapter Template

Every chapter `README.md` will use the following sections, omitting only a
section that genuinely does not apply:

1. **Overview** — the scientific or mathematical question and why it matters.
2. **Core ideas** — a plain-language explanation of the key mechanisms and
   terminology.
3. **Essential model** — the minimum equations, variables, parameters, and
   units needed to read the programs. Symbols must be defined immediately.
4. **Code examples** — one linked entry for every immediate example
   subdirectory, explaining what its program computes or plots.
5. **What to look for** — the important qualitative features in the generated
   traces, phase planes, bifurcation diagrams, raster plots, or other figures.
6. **Suggested order** — an ordering when the examples build on one another;
   otherwise, state that they can be read independently.
7. **Prerequisites and related chapters** — links to local chapter guides that
   introduce required or closely related concepts.
8. **Running the examples** — the normal command and any chapter-specific
   inputs, generated files, or runtime caveats.
9. **Further reading** — optional links or citations to suitable general and
   primary references.

The practical target is approximately one to two rendered Markdown pages per
chapter. Chapters with many examples may be longer because the code-example
index must remain complete, but conceptual exposition should stay concise.

## Content Method

For each chapter, inspect every immediate example directory and its main entry
point. Derive the example descriptions, variable names, parameter conventions,
and output interpretation from the repository implementation. Use the book only
to confirm the chapter's conceptual framing and sequence. Consult external
primary or authoritative references only when a concept cannot be explained
accurately from the implementation and supplied source material.

Descriptions must distinguish between:

- simulations of time-dependent neuronal or network models;
- analytical or numerical dynamical-systems calculations;
- parameter sweeps and derived curves;
- static explanatory plots.

Equations should be included only when they directly help a reader understand
the code. Each equation must use the same sign and variable conventions as the
implementation, or explicitly call out a local naming difference.

## Navigation and Running Instructions

Links must be relative so they work both in a local checkout and on GitHub.
Each chapter guide links to all immediate example directories, and the Python
index links to every chapter guide. Cross-chapter links should be limited to
genuine prerequisites or closely related follow-up material.

Running instructions will use the repository's current execution pattern:

```bash
cd python/<chapter>/<example>
python main.py
```

Individual guides will identify exceptions such as auxiliary plotting scripts,
input data, generated text files, longer network runs, or legacy makefiles.
Installation details remain centralized in the root documentation rather than
being copied into every chapter.

## Validation

Documentation validation will include:

- all 38 existing Python chapter directories contain a `README.md`;
- `python/README.md` links to all 38 guides and no nonexistent chapter;
- every immediate example subdirectory in a chapter is represented in that
  chapter's **Code examples** section;
- all relative links resolve;
- equations, parameter descriptions, commands, inputs, and outputs agree with
  the checked-in code;
- terminology and chapter mapping are consistent across the index and guides;
- no copied book passages, exercises, or book illustrations are included;
- a representative review of early, middle, and late chapters confirms the
  template is useful for single-neuron, dynamical-systems, network, and
  plasticity material.

Existing numerical tests need not change because this work does not alter
program behavior. A small documentation-check script or test may be added if it
is the simplest reliable way to enforce directory coverage and relative links.

## Scope Boundaries

- Document only the existing Python chapter directories and their current
  examples.
- Do not modify simulation behavior, regenerate figures, or reorganize example
  directories as part of this work.
- Do not add full textbook-style derivations, exercise sets, or broad
  neuroscience surveys.
- Do not duplicate a full guide under `brian/` or `matlab/`; those trees may
  link to the Python conceptual guides in later work.
- Preserve and do not stage unrelated uncommitted changes already present in
  the worktree.

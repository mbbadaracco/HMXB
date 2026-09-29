# Working agreement
Everything below is written for the agent. 

---

## The project

The person you are working with is an astronomer and does not need concepts explained unless she asks. The work we are about to begin consists on a project that must follow the indications in final-project.html (read and learn those rules). The first author is me, M. B. Badaracco, and my affiliation are:
1. Universidad de Buenos Aires, Facultad de Ciencias Exactas y Naturales, Departamento de Física. Buenos Aires, Argentina.
2. CONICET - Universidad de Buenos Aires, Instituto de Astronomía y Física del Espacio (IAFE). Buenos Aires, Argentina.
The context and goals of the project are the following. Read the papers in parenthesis, they are stored at /Papers:
Recent catalogs of Galactic HMXBs, Fortin2023 and Neumann2023, does not contain X-ray properties, besides a few like spin period. 
The Chandra Source Catalog 2.1.1 (Evans2024) as well as 5XMM-DR15 (5XMMdraft) does report X-ray properties of the HMXBs among so many other sources. However (CHECK THIS only by **reading** the papers of the telescope catalogs, not running anything) these works does not deconvolve the intrinsic absorption of the source from the absorption of the ISM. Both of these telescopes detect the same energy range, yet their catalogs report some of their characteristics in distinct wavebands, which makes them incomparable without further homogenization.
The main goal of the present work is to study the importance of incorporating absorption correction in the values derived from X-ray emission of HMXBs, through a new reduction of the data, for which I implement the 3D NH-tool developed by Doroshenko 2024 (Doroshenko2024), and a later comparison with the values reported by Chandra and XMM-Newton catalogs. Particular emphasis is put on homogenization between the observations of the two telescopes.

## Layout

```
5XMM-DR15/      latest version of the XMM-Newton catalog
Extra/          worry later about this
Papers/         source PDFs. Never cite a paper that is not in here
work/           scripts, notebooks, figures. All what YOU do lives here. Do NOT create or modify files outside this folder.
```

## How to work

- **Prefer a check to a claim.** If you derive something, verify it with `sympy` and assert the
  difference is zero. If you compute a number a paper printed, assert you reproduced it. A
  result with no assertion is a hypothesis; say so in those words.
- **Say which parts you derived and which you recalled.** Always, unprompted. They look
  identical on the page and only one of them survives a new problem.
- **Never claim to have run something you did not run.** "The script produces X" and "I ran the
  script and it printed X" are different sentences. Paste the output.
- **Reproducing a figure means reproducing the figure** — the same axes, the same curve, from
  the same kind of input — not describing it, and not drawing something with the right shape.
  If you could not get the data, say the figure was not reproduced.
- **State what you did *not* check.** A summary line like "12/12 passing" is misleading if half
  the work has no assertions. Name the uncovered half.

## Formats

- Write-ups: **LaTeX** compiled to PDF. Not markdown. If the software you need to compile latex is not installed, install it without asking.
- Figures: a **Jupyter** notebook (`.ipynb`), exported to HTML. 

## Do not

- Do not write long prose summaries. They are the one output that cannot be checked at a glance.
- Do not install anything into the base environment; make a per-exercise one.
- Do not create files outside work/.
- Do not `git commit` or `git push` unless asked.

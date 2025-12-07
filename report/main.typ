#import "@preview/elsearticle:1.1.0": *

/*
# The Research Catechism
1. Problem
 What is the problem? What opportunities does addressing it create?
2. Insight
 What is the key insight? What conceptual contributions does our solution provide?
3. Claims
 How do we know the key insight is good? What measurable claims can we make about the solution?
4. Evidence
 What evidence do we have to back up the claims?
*/

#let abstract = [
  this is a test of our abstract

  source code at https://github.com/srzer/cse-507-project
]

#show: elsearticle.with(
  title: "test",
  authors: (
    // a lipson, Jaedon Rich, Ruizhe Shi, Evan Wang
    (
      name: "Ruizhe Shi",
      affiliation: "University of Washington, Seattle, U.S.",
    ),
    (name: "Evan Wang"),
    (name: "Jaedon Rich"),
    (name: "lipson"),
  ),
  // journal: "UW CSE 507",
  abstract: abstract,
  keywords: ("interval arithmetic", "SMT"),
  format: "3p",
  // line-numbering: true,
)


/*
The purpose of the final report is to present your
- ideas,
- algorithms,
- implementation,
- and experiences—what worked and what didn’t.
It should again follow the narrative arc outlined in the research catechism and include the following sections:

The report should be 5 pages long (excluding references).
- Keep your writing brief and precise.
- Illustrate algorithms and encodings with examples and figures.
- Illustrate results with tables, graphs, and screenshots.
- Submit the report in PDF format by the due date.
- The final submission should also include
    - your demo slides and a
    - ZIP archive with the source code for the tool,
    - build scripts,
    - benchmarks (if applicable), and a
    - README file describing how to run the tool on your end-to-end scenarios.
*/

//  1. Introduction and precise problem statement.
= Introduction

#include "sections/introduction.typ"

// 2. Overview of your approach and a summary of how it relates to previous work.
= Overview

#include "sections/overview.typ"

//  3. Algorithms and Logical Encodings you developed.
= Algorithms

#include "sections/algorithms.typ"

/* 4. Summary of Results, including
key design and implementation challenges;
how you addressed them (what worked, what didn’t, and why);
and how this work could lead to a real tool or a full-length conference paper. */
= Results

#include "sections/summary.typ"

//  5. Teamwork: a one-paragraph description of the individual team member’s contributions.
= Teamwork

#include "sections/teamwork.typ"

/*  6. Course Topics: a one-paragraph description of the course topics applied in the project,
and a one-paragraph description (if applicable) of any topics that would have been useful but weren’t covered in the course. */
= Course Topics

#include "sections/course-topics.typ"

// = References (automatically included in template)
#bibliography("refs.bib")




/* template reference

#lorem(100)

= Section 1

#lorem(50)

== Subsection 1

#lorem(10) (see Eq. @eq1) @Aut10.

$
  y = alpha x + beta tau integral_0^x d x
$ <eq1>
where ...

$
  x = integral_0^x d x #<eqa>\
  (u v)' = u' v + v' u #<eqb>
$ <eq2>

Eq. @eqa is a simple integral, while Eq. @eqb is the derivative of a product of two functions. These equations are grouped in Eq. @eq2.

== Features

=== Table

Below is Table @tab:tab1.

#let tab1 = {
  table(
    columns: 3,
    table.header([*Header 1*], [*Header 2*], [*Header 3*]),
    [Row 1], [12.0], [92.1],
    [Row 2], [16.6], [104],
  )
}

#figure(tab1, kind: table, caption: [Example]) <tab:tab1>

=== Figures

Below is Fig. @fig:logo.

#figure(
  image("images/typst-logo.svg", width: 50%),
  caption: [Typst logo - Credit: \@fenjalien],
) <fig:logo>

=== Subfigures

Below are Figs. @figa and @figb, which are part of Fig. @fig:typst.

#subfigure(
  figure(image("images/typst-logo.svg"), caption: []),
  <figa>,
  figure(image("images/typst-logo.svg"), caption: []),
  <figb>,
  columns: (1fr, 1fr),
  caption: [(a) Left image and (b) Right image],
  label: <fig:typst>,
)

#show: appendix

= Appendix A

== Figures

In @fig:app

#figure(image("images/typst-logo.svg", width: 50%), caption: [Books cover]) <fig:app>

== Subfigures

Below are @figa-app and @figb-app, which are part of @fig:typst-app.

#subfigure(
  figure(image("images/typst-logo.svg"), caption: []),
  <figa-app>,
  figure(image("images/typst-logo.svg"), caption: []),
  <figb-app>,
  columns: (1fr, 1fr),
  caption: [(a) Left image and (b) Right image],
  label: <fig:typst-app>,
)

== Tables

In @tab:app

#figure(tab1, kind: table, caption: [Example]) <tab:app>

== Equations

In @eq

$
  y = f(x)
$ <eq>

#nonumeq[$
    y = g(x)
  $
]

$
  y = f(x) \
  y = g(x)
$

*/

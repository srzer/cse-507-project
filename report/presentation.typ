/*
Project Demo Evaluations

Mainly we want everyone to learn and be amazed by what folks have put together in teams over the quarter.
That is a bit vague though, so more concretely,
Audrey and I will be looking for the following in your demos for next Tuesday:

- Your tool works on a couple "engaging" end-to-end scenarios
- You give sufficient background for the problem and setting
- You clearly motivate why your project matters / is interesting
- You explain what your actual contribution is (not just what it does)
- Everyone on the team presents

We're really excited to see what you’ve built,
and the presentations are one of the best parts of the course.
Please show up ready to crush it and support your classmates' work too!
*/

#import "@preview/touying:0.6.1": *
#import "@preview/cetz:0.3.4" as cetz: *

#import themes.simple: *

#let cetz-canvas = touying-reducer.with(
  reduce: cetz.canvas,
  cover: cetz.draw.hide.with(bounds: true),
)

#show: simple-theme.with(aspect-ratio: "16-9")


= Hybrind Numerical & SMT Opmitization

Ruizhe Shi, Evan Wang, Jaedon Rich, a lipson

// The Research Catechism

== Problem
// What is the problem? What opportunities does addressing it create?

=== test

#figure(
  image("renders/examples_2.png", width: 50%),
  caption: [From left to right: hand-written example, Singular Edge, and Rational Bowl optimization search space renders.
  ],
)

== Insight
//What is the key insight? What conceptual contributions does our solution provide?

== Claims
// How do we know the key insight is good? What measurable claims can we make about the solution?

== Evidence
// What evidence do we have to back up the claims?/


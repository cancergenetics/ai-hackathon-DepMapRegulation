# M7 retro (fresh agent, process only; saved by the PI from the agent's reply)

**Stop doing**
1. Deferring LESSONS updates: none was written for D6–D8a. Write the line in the same commit as the decision note.
2. Spending tuning budget on objectives known to be confounded. Run cheap C1/C2 controls before adding a fusion tier, and don't write rules a variant cannot satisfy (I2).
3. Writing headline claims before verifying them. REPORT.md had 4 errors and 3 over-claims, and the research tree still carried them.
4. Patching rules late. D6a, the NaN CMH bound and the sealed preflight were caught by review or mock, not by tests. Require a known-answer toy test for each new metric.
5. Putting the Fisher p (2e-36) next to the pre-registered permutation p (0.068) in summaries.

**Worked well**
1. The final-run protocol: selection commit → LOCK commit → one completed run; TEST sealed since 57f2de3.
2. The pre-final fresh review and mock run: 3 real fixes, applied in about 13 minutes.
3. Evaluation-only controls (C1, C2, B0, NC1, NC2, class iii, B5) stopped over-interpretation of I1 (D7).
4. Independent number verification of REPORT.md (about 120 numbers checked).
5. Explicit disclosure of post-hoc work (D6), with lineage halves committed before evaluation.

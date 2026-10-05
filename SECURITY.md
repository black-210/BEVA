# Security

BEVA is a numerical library, not a cryptographic primitive.

For untrusted inputs:
- enforce application-level dimension limits;
- enforce coefficient limits;
- inspect `BEVA_ERR_OVERFLOW`;
- keep LLL `max_iterations` finite;
- do not treat floating-point LLL as an exact-arithmetic proof.

Report memory-safety or correctness issues through the project issue tracker.

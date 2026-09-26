#!/usr/bin/env python3
"""Conservative explicit perturbation windows for square Ising blocks.

All input temperatures and certified margins are exact decimal strings.
Positive upper bounds use decimal ROUND_CEILING; lower bounds use
ROUND_FLOOR.  Decimal exp/ln are correctly rounded; one adjacent decimal
is added to enclose their exact values.  No binary floating point enters
the certificate.  This script does not certify the supplied Ising margin.

The shell bound is C_s (r+1)^3 exp(-mu*r), with the constants derived in
the accompanying certification_details.tex shell calculation.  A geometric
majorant controls the infinite tail, rather than truncating that tail.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext


D = Decimal
PRECISION = 100


def lower_div(a: Decimal, b: Decimal) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_FLOOR
        return a / b


def lower_mul(a: Decimal, b: Decimal) -> Decimal:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_FLOOR
        return a * b


def upper_exp(a: Decimal) -> Decimal:
    # Decimal.exp is correctly rounded with ROUND_HALF_EVEN irrespective
    # of the context rounding mode.  next_plus encloses that last ulp.
    return a.exp().next_plus()


def upper_ln(a: Decimal) -> Decimal:
    return a.ln().next_plus()


def window(x_string: str, ell: int, kappa_string: str, u0_string: str) -> dict:
    with localcontext() as ctx:
        ctx.prec = PRECISION
        ctx.rounding = ROUND_CEILING
        x, kappa, u0 = map(D, (x_string, kappa_string, u0_string))
        if not (x > 0 and ell >= 1 and 0 < kappa <= 1 and 0 < u0 <= 1):
            raise ValueError("Require x>0, ell>=1, 0<kappa<=1, 0<u0<=1")

        b = D(ell * ell)
        r0 = D(ell * ell + 4 * ell)
        chi = r0 / b
        volume, density = r0, chi
        w = D(4 * ell) * x
        E = upper_exp(w)
        a0 = r0 / 2
        projector_count = D(2) ** int(r0)
        pin_bound = projector_count * a0 * (2 + a0 * u0)
        A_T = E * (3 + 2 * w)
        C_loc = 12 * A_T + E / 2 + 12 * pin_bound
        T = E * upper_ln(32 * r0 * chi / kappa)

        C_LR = D(8)
        vbar = max(D(1), 32 * upper_exp(D(1)) * x)
        mu = lower_div(D(1), 2 * vbar)
        C_transform = 4 * C_LR + 4 * vbar + 4
        K_Q = E * C_transform * r0 * w
        L_Q = max(E, upper_exp(mu) * K_Q)
        jump_bound = 1 + a0 * u0
        K_P = 4 * projector_count * jump_bound * C_transform * r0
        L_P = max(pin_bound * u0, upper_exp(mu) * K_P)
        C_s = ((8 * r0 + 34) * (K_Q + K_P)
               + 2 * C_transform * r0 * (L_Q + L_P))
        A = 2 * T * C_loc * volume ** 3 * density
        B0 = 2 * T * C_s * volume * density

        # For R+2>=14/mu, t_(r+1)/t_r <= exp(-mu/2), where
        # t_r=(r+2)^7 exp(-mu*r).  Therefore the whole tail is at most
        # (R+2)^7 exp(-mu*R)/(1-exp(-mu/2)).
        half_mu_lower = lower_div(mu, D(2))
        q_upper = upper_exp(-half_mu_lower)
        with localcontext() as down:
            down.prec = PRECISION
            down.rounding = ROUND_FLOOR
            tail_denominator_lower = D(1) - q_upper
            target_lower = kappa / 8
        if tail_denominator_lower <= 0:
            raise ArithmeticError("Insufficient precision for tail denominator")

        def tail_upper(radius: int) -> Decimal:
            decay = upper_exp(-lower_mul(mu, D(radius)))
            return (B0 * D(radius + 2) ** 7 * decay
                    / tail_denominator_lower)

        radius_min = max(0, int((14 / mu - 2).to_integral_value(
            rounding=ROUND_CEILING)))
        radius_hi = max(1, radius_min)
        while tail_upper(radius_hi) > target_lower:
            radius_hi *= 2
        radius_lo = radius_min - 1
        while radius_hi - radius_lo > 1:
            mid = (radius_hi + radius_lo) // 2
            if tail_upper(mid) <= target_lower:
                radius_hi = mid
            else:
                radius_lo = mid
        radius = radius_hi
        u_lower = min(u0, lower_div(kappa, 8 * A * D(1 + radius) ** 8))
        # A power of ten strictly inside the sufficient interval is easy
        # to quote and leaves ample slack for a reader's reproduction.
        u_safe = D(1).scaleb(u_lower.adjusted() - 1)
        local_error_upper = A * u_safe * D(1 + radius) ** 8
        certified_tail_upper = tail_upper(radius)
        assert local_error_upper <= target_lower
        assert certified_tail_upper <= target_lower

        values = {
            "x": x, "ell": ell, "kappa_lower_input": kappa, "u0": u0,
            "b": b, "r0": r0, "chi_upper": chi,
            "w": w, "E_upper": E, "T_upper": T,
            "C_loc_upper": C_loc, "C_s_upper": C_s,
            "mu_lower": mu, "A_upper": A, "B0_upper": B0,
            "R_star": radius, "tail_error_upper": certified_tail_upper,
            "local_error_upper_at_safe_u": local_error_upper,
            "u_star_lower": u_lower, "safe_beta_abs_g": u_safe,
        }
        return {key: str(value) if isinstance(value, D) else value
                for key, value in values.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x", help="Exact decimal beta*J")
    parser.add_argument("--ell", type=int, help="Square side length")
    parser.add_argument("--kappa", help="Certified lower bound on margin")
    parser.add_argument("--u0", default="1", help="Initial perturbation cap")
    args = parser.parse_args()
    if any(a is not None for a in (args.x, args.ell, args.kappa)):
        if any(a is None for a in (args.x, args.ell, args.kappa)):
            parser.error("--x, --ell, --kappa must be supplied together")
        cases = [(args.x, args.ell, args.kappa)]
    else:
        cases = [("0.35", 6, "0.0432"), ("0.36", 7, "0.0125")]
    print(json.dumps([window(x, ell, kappa, args.u0)
                      for x, ell, kappa in cases], indent=2))


if __name__ == "__main__":
    main()

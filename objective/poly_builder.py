from collections import defaultdict
from .polynomial_type import Polynomial, to_term
from .rational_type import Rational


class PolyBuilder:
    """
    Helper class to build Polynomials and Rationals using natural python syntax
    (e.g., x**2 + y).
    """

    def __init__(self, terms=None, dim=3):
        self.dim = dim
        self.terms = defaultdict(float)
        if terms:
            for exps, coeff in terms.items():
                self.terms[exps] = coeff

    @classmethod
    def var(cls, idx, dim=3):
        exps = [0] * dim
        exps[idx] = 1
        return cls({tuple(exps): 1.0}, dim)

    @classmethod
    def const(cls, val, dim=3):
        exps = [0] * dim
        return cls({tuple(exps): float(val)}, dim)

    def __add__(self, other):
        other = self._ensure_poly(other)
        new_terms = self.terms.copy()
        for exps, coeff in other.terms.items():
            new_terms[exps] += coeff
        return PolyBuilder(new_terms, self.dim)

    def __sub__(self, other):
        other = self._ensure_poly(other)
        new_terms = self.terms.copy()
        for exps, coeff in other.terms.items():
            new_terms[exps] -= coeff
        return PolyBuilder(new_terms, self.dim)

    def __mul__(self, other):
        other = self._ensure_poly(other)
        new_terms = defaultdict(float)
        for e1, c1 in self.terms.items():
            for e2, c2 in other.terms.items():
                new_exps = tuple(a + b for a, b in zip(e1, e2))
                new_terms[new_exps] += c1 * c2
        return PolyBuilder(new_terms, self.dim)

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("Only non-negative integer powers supported")
        res = PolyBuilder.const(1.0, self.dim)
        base = self
        for _ in range(power):
            res = res * base
        return res

    def __truediv__(self, other):
        other = self._ensure_poly(other)
        return Rational(self.to_polynomial(), other.to_polynomial())

    def _ensure_poly(self, other):
        if isinstance(other, (int, float)):
            return PolyBuilder.const(other, self.dim)
        return other

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return PolyBuilder.const(other, self.dim) - self

    def __rmul__(self, other):
        return self * other

    def to_polynomial(self) -> Polynomial:
        term_list = []
        for exps, coeff in self.terms.items():
            if abs(coeff) > 1e-9:
                term_list.append(to_term(coeff, list(exps)))
        # If empty (zero), return a zero term
        if not term_list:
            term_list.append(to_term(0.0, [0] * self.dim))
        return Polynomial(term_list)


def make_vars(dim=3):
    return [PolyBuilder.var(i, dim) for i in range(dim)]

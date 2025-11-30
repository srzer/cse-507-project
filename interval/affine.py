from dataclasses import dataclass


@dataclass(frozen=True)
class AffineForm:
    """
    Very lightweight 'affine' form:
      x ≈ mid + sum_i coeffs[i] * ε_i,  with ε_i ∈ [-1,1].
    We are not doing full affine propagation here; we only use
    mid/coeffs to get a cheap range and rescale for sub-boxes.
    """

    mid: float
    coeffs: list[float]

    def range(self):
        r = sum(abs(c) for c in self.coeffs)
        return (self.mid - r, self.mid + r)

    def scale_shift(self, a_new, b_new):
        """
        Rescale this affine form from its current interval to [a_new, b_new].
        This is used when we subdivide a box and want to reuse parent AA info.
        """
        old_min, old_max = self.range()
        r_width = old_max - old_min
        r_new_width = b_new - a_new
        if r_width == 0:
            # degenerate: collapse to midpoint of new interval
            alpha = 0.0
            beta = (a_new + b_new) * 0.5
        else:
            alpha = r_new_width / r_width
            beta = a_new - alpha * old_min

        new_coeffs = [c * alpha for c in self.coeffs]
        new_mid = self.mid * alpha + beta
        return AffineForm(new_mid, new_coeffs)

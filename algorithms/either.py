from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Any, Callable, Generic, TypeVar
from typing_extensions import Self

V = TypeVar("V")
E = TypeVar("E")
U = TypeVar("U")


@dataclass(frozen=True)
class Left(Generic[E]):
    error: E

    # just pass the error along; f is never called
    def map(self: Self, f: Callable[[Any], Any]) -> Self:
        return self

    def bind(self: Self, f: Callable[[Any], Any]) -> Self:
        return self

    def fold(self, l: Callable[[E], U], r: Callable[[Any], U]) -> U:
        return l(self.error)

    def get_or_else(self, default: U) -> U:
        return default


@dataclass(frozen=True)
class Right(Generic[V]):
    value: V

    def map(self, f: Callable[[V], U]) -> Self:
        return replace(self, value=Right(f(self.value)))

    def bind(self, f: Callable[[V], Either]) -> Either:
        return f(self.value)

    def fold(self, l: Callable[[Any], U], r: Callable[[V], U]) -> U:
        return r(self.value)

    def get_or_else(self, default: Any) -> V:
        return self.value


Either = Left[E] | Right[V]


(** This file is part of RAGlib (Real Algebraic Geometry Library).
 *
 * RAGlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * RAGlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RAGlib.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Mohab Safey El Din **)

#Exported functions of RAGlib

HasRealSolutions:=proc(eqs, pos, ineqs, opts:={})
local isempty, _l, newopts;

  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty");
  else 
    isempty:=0:
  end if;

  if isempty = 0 then 
    return SemiAlgebraicSolve(eqs, pos, ineqs, opts union {"isempty"=1});
  else 
    newopts:=map(_l->if lhs(_l)="isempty" then lhs(_l)=1 else _l fi, opts);
    return SemiAlgebraicSolve(eqs, pos, ineqs, newopts);
  end if;
end proc:

PointsPerComponents:=proc(eqs, pos, ineqs, opts:={})
local isempty, _l, newopts;

  if type(subs(opts, "isempty"), integer) then 
    isempty:=subs(opts, "isempty");
  else 
    isempty:=0:
  end if;

  if isempty = 0 then 
    return SemiAlgebraicSolve(eqs, pos, ineqs, opts);
  else 
    newopts:=map(_l->if lhs(_l)="isempty" then lhs(_l)=0 else _l fi, opts);
    return SemiAlgebraicSolve(eqs, pos, ineqs, newopts);
  end if;
end proc:

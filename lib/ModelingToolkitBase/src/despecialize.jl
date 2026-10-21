struct __Despecialize1 <: ReversibleTransformations end
struct __Despecialize2 <: ReversibleTransformations end
struct __Despecialize3 <: ReversibleTransformations end

reverse_transformation_by_default(::__Despecialize1) = true
reverse_transformation_by_default(::__Despecialize2) = true
reverse_transformation_by_default(::__Despecialize3) = true

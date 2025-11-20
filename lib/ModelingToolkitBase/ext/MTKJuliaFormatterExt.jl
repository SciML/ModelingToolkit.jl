module MTKJuliaFormatterExt

import ModelingToolkitBase: readable_code, _readable_code, rec_remove_macro_linenums!
import JuliaFormatter

function readable_code(expr::Expr)
    expr = Base.remove_linenums!(_readable_code(expr))
    rec_remove_macro_linenums!(expr)
    JuliaFormatter.format_text(string(expr), JuliaFormatter.SciMLStyle())
end

end

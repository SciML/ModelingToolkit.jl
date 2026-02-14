SimpleArrayModel.fmu is a v2 ME FMU that contains multidimensional array variables with derivatives. It is included to test SciML issue number 3934. It is generated from the Modelica model in test_mtk_2darray_bug.mo using  Dymola 2022x. The translation command is

```julia
translateModelFMU("SimpleArrayModel", false, "", "2", "me", true, 0, fill("",0));
```
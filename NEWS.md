# ModelingToolkit v8 Release Notes


### Upgrade guide

- `connect` should not be overloaded by users anymore. `[connect = Flow]`
  informs MTK that currents ought to sum to zero, and by default, variables are
  equal in a connection. Please check out [acausal components tutorial](https://mtk.sciml.ai/dev/tutorials/acausal_components/)
  for examples.

# ModelingToolkit v8 Release Notes


### Upgrade guide

- `connect` should not be overloaded by users anymore. `[connect = Flow]`
  informs ModelingToolkit that particular variable in a connector ought to sum
  to zero, and by default, variables are equal in a connection. Please check out
  [acausal components tutorial](https://docs.sciml.ai/ModelingToolkit/stable/tutorials/acausal_components/)
  for examples.

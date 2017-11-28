cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement

inputs: []

outputs: []

baseCommand: [echo, -n]

arguments:
  - "Stepik course"
  - "echo"
  - |
    CWL is not so simple )
  - "|"
  - "wc"
  - "-c"

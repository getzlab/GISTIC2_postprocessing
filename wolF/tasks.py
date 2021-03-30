from wolf import Task

def task1(
  <task_arg1> = <arg1_default>,
  <task_arg2> = <arg2_default>,
  ...,
  **task_kwargs
):
    # preprocessing logic can go here

    return Task(
      name = <task name>,
      inputs = {
        "arg1" : <task_arg1>,
        "arg2" : <task_arg2>,
        ...
      },
      script = [
        """
        # bash script goes here!
        # $arg1 and $arg2 are available as variables to use, e.g.
        echo $arg1 $arg2
        """
      ],
      outputs = {
        "output1" : "<output1 pattern>",
        "output2" : "<output2 pattern>",
        ...
      },
      docker = "gcr.io/broad-getzlab-workflows/<container name>:<container version>",
      resources = { "mem" : "1G" },
      **task_kwargs
    )

# define additional tasks in the same way that task1 is defined above.

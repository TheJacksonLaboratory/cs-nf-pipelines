process RENAME{

  input:
  tuple val(sampleID), file(file)
  val(new_name)

  output:
  tuple val(sampleID), file("${new_name}")

  script:
  """
  mv ${file} ${new_name}
  """

}

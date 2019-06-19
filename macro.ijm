function action(input, output, filename) {
  open(input + filename);
  run("Subtract Background...","rolling=50");
  //run("Subtract Background...", "rolling=50 stack");
saveAs("Tiff", output + filename);
close();
  print(input + filename);
}

#@String input
#@String output
extension2 = ".tiff";


list = getFileList(input);
for (i=0; i < list.length; i++)
  action(input, output, list[i]);

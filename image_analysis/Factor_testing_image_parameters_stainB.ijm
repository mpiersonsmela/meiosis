n = nImages;
print(n);
for (i=1;i<=nImages;i++)
{
print(i);
selectImage(i);
originalName = getTitle();
print(originalName);
run("Scale Bar...", "width=50 height=40 thickness=20 font=100 color=White background=None location=[Lower Right] horizontal bold overlay");
//run("Brightness/Contrast...");
setMinAndMax(256, 46829);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_DAPI_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(0, 65279);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_SYCP3_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 65279);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_Phalloidin_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(768, 16720);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_TEX12_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(65535, 65535);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_DMC1_.png");
run("Make Composite");
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_stainB/F3/"+originalName+"_overlay_.png");
}

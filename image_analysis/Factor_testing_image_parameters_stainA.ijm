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
setMinAndMax(256, 40807);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_DAPI_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(768, 31454);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_SYCP3_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 48768);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_Phalloidin_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(768, 44779);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_HORMAD1_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 60859);
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_gammaH2AX_.png");
run("Make Composite");
saveAs("PNG", "C:/Users/Wyss User/Pictures/TF_testing_images/F2/"+originalName+"_overlay_.png");
}
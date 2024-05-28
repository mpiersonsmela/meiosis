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
setMinAndMax(768, 48815);
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_HORMAD1_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 65535);
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_gammaH2X_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(0, 65535);
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_DAPI_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(768, 62204);
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_SYCP3_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 62844);
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_Phalloidin_.png");
run("Make Composite");
saveAs("PNG", "C:/Users/Wyss User/Pictures/2024_03_29_temp_testing/Slide4/"+originalName+"_overlay_.png");
}
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
setMinAndMax(256, 46124);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_DAPI_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(512, 52146);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_SYCP3_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 32768);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_T2A_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(0, 24087);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_HORMAD1_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 40807);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_gammaH2AX_.png");
run("Make Composite");
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainA/F2_2/"+originalName+"_overlay_.png");
}
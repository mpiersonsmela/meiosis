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
setMinAndMax(256, 65535);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_DAPI_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(0, 61499);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_SYCP3_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 54196);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_KI67_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(0, 26778);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_TEX12_.png");
run("Next Slice [>]");
//run("Brightness/Contrast...");
setMinAndMax(256, 65535);
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_RAD51_.png");
run("Make Composite");
saveAs("PNG", "C:/Users/Wyss User/Pictures/Meiosis_15_day_timecourse_stainB/F2/"+originalName+"_overlay_.png");
}

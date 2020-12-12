

#cp 84-J238-3-1_27F_TSS20201125-025-7163_D10.seq 84-J238-3-1_27F_TSS20201125-025-7163_D10.txt

#sed -i 's/\r//g' 84-J238-3-1_27F_TSS20201125-025-7163_D10.txt

#echo `tr -d "\n"  <84-J238-3-1_27F_TSS20201125-025-7163_D10.txt ` >84-J238-3-1_27F_TSS20201125-025-7163_D10.txt
#------------------------------------------------------------------



ls *.seq > count.list
sed 's/.seq//g' count.list > txt.list

for i in `cat txt.list` ;do
  cp ${i}.seq  ${i}.txt
  sed -i 's/\r//g' ${i}.txt
  echo `tr -d "\n"  <${i}.txt ` >${i}.txt
done

rm *.list





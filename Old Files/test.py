import matplotlib.pyplot as plt
import pandas as pd
  
df = pd.read_csv("POWERFACTOR.csv")
print(df)
# Extract the data from row 2046
row_data = df.iloc[2045:] # Rows are indexed from 0, so row 2046 is indexed 2045

# Set the x-axis data (Time, s)
x_data = row_data['Time, s']

# Set the y-axis data for the 2nd and 5th columns
y1_data = row_data['VPS1 [Case 1]']
y2_data = row_data['CPS1 [Case 1]']

# Create a plot
plt.figure(figsize=(10, 6))
plt.plot([x_data], [y1_data], marker='o', label='VPS1 [Case 1]', linestyle='None')
plt.plot([x_data], [y2_data], marker='o', label='CPS1 [Case 1]', linestyle='None')

# Add labels and a legend
plt.xlabel('Time, s')
plt.ylabel('Value')
plt.legend()

# Display the plot
plt.show()

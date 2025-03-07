FROM openjdk:8u102

# Download the tool
RUN wget https://github.com/dieterich-lab/JACUSA2/releases/download/v2.0.4/JACUSA_v2.0.4.jar -O /usr/local/bin/JACUSA_v2.0.4.jar

# Create a wrapper for easy access
RUN echo "#!/usr/bin/env bash" > /usr/local/bin/jacusa2
RUN echo "java -jar /usr/local/bin/JACUSA_v2.0.4.jar $@" >> /usr/local/bin/jacusa2
RUN chmod +x /usr/local/bin/jacusa2

FROM minio/minio:RELEASE.2022-01-07T01-53-23Z
COPY --from=minio/mc:RELEASE.2022-01-07T06-01-38Z /bin/mc /opt/bin/mc

EXPOSE 9000
EXPOSE 9001

ADD ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT [ "/bin/sh", "-c" ]
CMD [ "/entrypoint.sh" ]

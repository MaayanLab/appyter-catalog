FROM minio/minio:RELEASE.2022-04-12T06-55-35Z
COPY --from=minio/mc:RELEASE.2022-04-07T21-43-27Z /bin/mc /opt/bin/mc

EXPOSE 9000
EXPOSE 9001

ADD ./entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT [ "/bin/sh", "-c" ]
CMD [ "/entrypoint.sh" ]
